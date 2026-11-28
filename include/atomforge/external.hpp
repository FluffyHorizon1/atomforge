// atomforge — ExternalPotential: energies, forces and virials from a persistent worker process
//
// The worker is spawned once (through /bin/sh -c "exec <command>"; without "exec" when the command uses shell operators or starts with a builtin) and talks over its stdin/stdout, one JSON object per line.
// Its stderr is appended to a log file. The full protocol is documented in docs/EXTERNAL_PROTOCOL.md; in short:
//
//   -> {"op":"hello","protocol":1,"elements":["Al"]}
//   <- {"protocol":1,"name":"pyscf(PBE/gth-szv)","cutoff":0,"masses":[26.98]}        ("masses" and "cutoff" optional)
//   -> {"op":"evaluate","id":1,"natoms":N,"cell":[[..],[..],[..]],"positions":[[x,y,z],...],"species":[0,...],"symbols":["Al",...]}
//   <- {"id":1,"energy":E,"forces":[[fx,fy,fz],...],"virial":[[..],[..],[..]]}       (or "stress" instead of "virial")
//   <- {"id":1,"error":"message"}                                                     (worker-side failure)
//   -> {"op":"exit"}
//
// Units: A, eV, eV/A. "virial" is W_ab = dE/d(eps_ab) in eV, the convention of ZhouEAM (tensile positive: W > 0 when
// stretching raises the energy). "stress" is W / V in eV/A^3 and is converted with the current cell volume.
// Cell rows are lattice vectors, as in Structure.
//
// Failure handling: a worker that cannot be started, exits, closes its stdout, exceeds the timeout, answers with malformed
// JSON, a wrong id, missing fields, the wrong number of atoms or vectors, or any non-finite number makes evaluate() throw
// ExternalPotentialError. The message includes the tail of the worker's stderr log. A worker that failed is killed and is not
// restarted: later calls throw too. A reply {"id":N,"error":"..."} also throws, but leaves the worker running (an SCF that
// did not converge should not end a search). POSIX systems with /bin/sh (Linux, macOS; tested on Linux): fork/exec,
// socketpair, pipe, poll. The worker runs under an exec'd /bin/sh supervisor (see supervisor_script()).
// One request at a time: a concurrent evaluate() from a second thread throws.
#pragma once
#include "potential.hpp"
#include <atomic>
#include <cerrno>
#include <charconv>
#include <locale>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <memory>
#include <poll.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <sys/socket.h>
#if defined(__linux__)
#include <sys/syscall.h>
#endif
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#if defined(MSG_NOSIGNAL)
#define ATOMFORGE_MSG_NOSIGNAL MSG_NOSIGNAL
#else
#define ATOMFORGE_MSG_NOSIGNAL 0   // e.g. macOS: SO_NOSIGPIPE is set on the socket instead
#endif
// Floating-point std::from_chars / std::to_chars (libstdc++ 11+, MSVC, recent libc++); otherwise streams in the classic
// locale, which are just as locale independent.
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
#define ATOMFORGE_FP_CHARCONV 1
#else
#define ATOMFORGE_FP_CHARCONV 0
#endif

namespace atomforge {

inline constexpr int external_protocol_version = 1;

struct ExternalPotentialError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// ---------------------------------------------------------------- minimal JSON (enough for the protocol)
namespace json {

struct Value {
    enum class Type { Null, Bool, Number, String, Array, Object } type = Type::Null;
    bool b = false;
    double num = 0;
    std::string str;
    std::vector<Value> arr;
    std::map<std::string, Value> obj;

    bool is_number() const { return type == Type::Number; }
    bool is_array() const { return type == Type::Array; }
    bool is_object() const { return type == Type::Object; }
    bool is_string() const { return type == Type::String; }
    const Value* get(const std::string& k) const {
        if (type != Type::Object) return nullptr;
        auto it = obj.find(k);
        return it == obj.end() ? nullptr : &it->second;
    }
};

class Parser {
public:
    explicit Parser(const std::string& s) : s_(s) {}
    Value parse() {
        Value v = value();
        ws();
        if (i_ != s_.size()) fail("trailing characters");
        return v;
    }

private:
    const std::string& s_;
    std::size_t i_ = 0;
    int depth_ = 0;

    [[noreturn]] void fail(const std::string& what) const {
        throw ExternalPotentialError("malformed JSON from worker (" + what + " at byte " + std::to_string(i_) + ")");
    }
    void ws() { while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t' || s_[i_] == '\r' || s_[i_] == '\n')) ++i_; }
    bool lit(const char* w) {
        const std::size_t n = std::strlen(w);
        if (s_.compare(i_, n, w) == 0) { i_ += n; return true; }
        return false;
    }
    Value value() {
        if (++depth_ > 64) fail("nesting too deep");
        ws();
        if (i_ >= s_.size()) fail("unexpected end");
        Value v;
        const char c = s_[i_];
        if (c == '{') {
            v.type = Value::Type::Object; ++i_; ws();
            if (i_ < s_.size() && s_[i_] == '}') { ++i_; --depth_; return v; }
            for (;;) {
                ws();
                if (i_ >= s_.size() || s_[i_] != '"') fail("expected key");
                std::string k = string();
                ws();
                if (i_ >= s_.size() || s_[i_] != ':') fail("expected ':'");
                ++i_;
                v.obj[k] = value();
                ws();
                if (i_ < s_.size() && s_[i_] == ',') { ++i_; continue; }
                if (i_ < s_.size() && s_[i_] == '}') { ++i_; break; }
                fail("expected ',' or '}'");
            }
        } else if (c == '[') {
            v.type = Value::Type::Array; ++i_; ws();
            if (i_ < s_.size() && s_[i_] == ']') { ++i_; --depth_; return v; }
            for (;;) {
                v.arr.push_back(value());
                ws();
                if (i_ < s_.size() && s_[i_] == ',') { ++i_; continue; }
                if (i_ < s_.size() && s_[i_] == ']') { ++i_; break; }
                fail("expected ',' or ']'");
            }
        } else if (c == '"') {
            v.type = Value::Type::String; v.str = string();
        } else if (lit("true")) { v.type = Value::Type::Bool; v.b = true; }
        else if (lit("false")) { v.type = Value::Type::Bool; v.b = false; }
        else if (lit("null")) { v.type = Value::Type::Null; }
        // Python's json module writes NaN / Infinity for non-finite floats. Accept them here so the caller can report
        // "non-finite value" rather than a parse error.
        else if (lit("NaN")) { v.type = Value::Type::Number; v.num = std::nan(""); }
        else if (lit("Infinity")) { v.type = Value::Type::Number; v.num = HUGE_VAL; }
        else if (lit("-Infinity")) { v.type = Value::Type::Number; v.num = -HUGE_VAL; }
        else if (c == '-' || (c >= '0' && c <= '9')) {
            // std::from_chars: independent of the C locale (a "de_DE" decimal comma must not break the protocol)
            const char* b = s_.data() + i_;
            const char* e = s_.data() + s_.size();
            v.type = Value::Type::Number;
#if ATOMFORGE_FP_CHARCONV
            const auto res = std::from_chars(b, e, v.num, std::chars_format::general);
            if (res.ec != std::errc() || res.ptr == b) fail("bad number '" + s_.substr(i_, 24) + "'");
            i_ += std::size_t(res.ptr - b);
#else
            std::size_t n = 0;
            while (b + n < e && std::strchr("+-.eE0123456789", b[n]) && b[n] != '\0') ++n;
            std::istringstream is(std::string(b, n));
            is.imbue(std::locale::classic());
            is >> v.num;
            if (!is || is.peek() != std::char_traits<char>::eof()) fail("bad number '" + s_.substr(i_, 24) + "'");
            i_ += n;
#endif
            if (i_ < s_.size() && s_[i_] != '\0' && std::strchr(".eE+-0123456789", s_[i_])) fail("bad number '" + s_.substr(std::size_t(b - s_.data()), 24) + "'");
        } else fail(std::string("unexpected character '") + c + "'");
        --depth_;
        return v;
    }
    std::string string() {
        ++i_;   // opening quote
        std::string out;
        while (i_ < s_.size() && s_[i_] != '"') {
            char c = s_[i_++];
            if (c == '\\') {
                if (i_ >= s_.size()) fail("bad escape");
                const char e = s_[i_++];
                switch (e) {
                    case '"': out += '"'; break;
                    case '\\': out += '\\'; break;
                    case '/': out += '/'; break;
                    case 'b': out += '\b'; break;
                    case 'f': out += '\f'; break;
                    case 'n': out += '\n'; break;
                    case 'r': out += '\r'; break;
                    case 't': out += '\t'; break;
                    case 'u': {   // keep ASCII, replace anything else by '?': only used in error messages
                        if (i_ + 4 > s_.size()) fail("bad \\u escape");
                        const unsigned long cp = std::strtoul(s_.substr(i_, 4).c_str(), nullptr, 16);
                        out += cp < 128 ? char(cp) : '?';
                        i_ += 4;
                        break;
                    }
                    default: fail("bad escape");
                }
            } else out += c;
        }
        if (i_ >= s_.size()) fail("unterminated string");
        ++i_;
        return out;
    }
};

inline Value parse(const std::string& s) { return Parser(s).parse(); }

inline std::string quote(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if (static_cast<unsigned char>(c) < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); o += b; }
        else o += c;
    }
    return o + "\"";
}
// Shortest representation that round-trips (std::to_chars: locale independent).
inline std::string num(double x) {
#if ATOMFORGE_FP_CHARCONV
    char b[64];
    const auto res = std::to_chars(b, b + sizeof b, x);
    return std::string(b, res.ptr);
#else
    std::ostringstream o;
    o.imbue(std::locale::classic());
    o.precision(17);
    o << x;
    return o.str();
#endif
}

} // namespace json

// ---------------------------------------------------------------- the potential
struct ExternalOptions {
    double hello_timeout_s = 60;     // start-up and handshake (imports of PySCF / torch can be slow)
    double eval_timeout_s = 3600;    // one evaluate request; DFT evaluations can take minutes
    std::string log_path;            // worker stderr; empty -> $TMPDIR/atomforge-worker-<pid>-<n>.log, one per instance
                                     // (timeouts >= 1e9 s, or infinite, mean no deadline)
    std::vector<double> masses;      // per element; overrides the worker's "masses" (else Zhou04 table, else error)
};

class ExternalPotential final : public Potential {
public:
    ExternalPotential(std::string command, std::vector<std::string> els, ExternalOptions opt = {})
        : cmd_(std::move(command)), names_(std::move(els)), opt_(std::move(opt)) {
        if (names_.empty()) throw ExternalPotentialError("ExternalPotential: no elements");
        if (opt_.log_path.empty()) opt_.log_path = default_log_path();
        start();
        handshake();
    }
    ~ExternalPotential() override { shutdown(); }
    ExternalPotential(const ExternalPotential&) = delete;
    ExternalPotential& operator=(const ExternalPotential&) = delete;

    double cutoff() const override { return cutoff_; }
    const std::vector<std::string>& elements() const override { return names_; }
    double mass(int s) const override {
        if (s < 0 || std::size_t(s) >= masses_.size() || !(masses_[s] > 0))
            throw ExternalPotentialError("ExternalPotential: no mass for species " + std::to_string(s) + " (pass ExternalOptions::masses)");
        return masses_[s];
    }
    std::string name() const override { return "external:" + worker_name_; }
    const std::string& log_path() const { return opt_.log_path; }
    long evaluations() const { return next_id_ - 1; }
    bool alive() const { return pid_ > 0; }
    pid_t supervisor_pid() const { return pid_; }   // leader of the worker's process group (for tests and diagnostics)

    EvalResult evaluate(const Structure& s, NeighborList&) const override {
        // One worker serves one request at a time. atomforge's relax / search / phonons call evaluate() sequentially; a
        // second, concurrent call (another thread) is refused rather than interleaved on the pipe.
        if (busy_.exchange(true)) throw ExternalPotentialError("ExternalPotential::evaluate called concurrently; use one ExternalPotential (worker) per thread");
        struct Release { std::atomic<bool>& b; ~Release() { b = false; } } release{busy_};
        if (pid_ <= 0) throw ExternalPotentialError("external worker is not running (it failed earlier; see " + opt_.log_path + ")");
        const std::size_t N = s.size();
        if (s.species.size() != N) throw ExternalPotentialError("structure has " + std::to_string(s.species.size()) + " species for " + std::to_string(N) + " atoms");
        for (int k = 0; k < 9; ++k)
            if (!std::isfinite(s.cell.m[k])) throw ExternalPotentialError("cell entry is not finite; not sent to the worker");
        for (std::size_t i = 0; i < N; ++i)
            if (!std::isfinite(s.pos[i].x) || !std::isfinite(s.pos[i].y) || !std::isfinite(s.pos[i].z))
                throw ExternalPotentialError("position of atom " + std::to_string(i) + " is not finite; not sent to the worker");
        const long id = next_id_++;
        std::string q;
        q.reserve(128 + 80 * N);
        q += "{\"op\":\"evaluate\",\"id\":" + std::to_string(id) + ",\"natoms\":" + std::to_string(N) + ",\"cell\":[";
        for (int i = 0; i < 3; ++i) q += std::string(i ? "," : "") + "[" + json::num(s.cell(i, 0)) + "," + json::num(s.cell(i, 1)) + "," + json::num(s.cell(i, 2)) + "]";
        q += "],\"positions\":[";
        for (std::size_t i = 0; i < N; ++i) q += std::string(i ? "," : "") + "[" + json::num(s.pos[i].x) + "," + json::num(s.pos[i].y) + "," + json::num(s.pos[i].z) + "]";
        q += "],\"species\":[";
        for (std::size_t i = 0; i < N; ++i) {
            const int sp = s.species[i];
            if (sp < 0 || std::size_t(sp) >= names_.size()) throw ExternalPotentialError("species index " + std::to_string(sp) + " out of range");
            q += (i ? "," : "") + std::to_string(sp);
        }
        q += "],\"symbols\":[";
        for (std::size_t i = 0; i < N; ++i) q += (i ? "," : "") + json::quote(names_[s.species[i]]);
        q += "]}\n";

        const json::Value r = request(q, opt_.eval_timeout_s, "evaluate");
        // A reported error (e.g. an SCF that did not converge) leaves the worker usable for the next structure.
        if (const json::Value* e = r.get("error"); e && e->type != json::Value::Type::Null && r.get("id") && r.get("id")->is_number() && r.get("id")->num == double(id))
            throw ExternalPotentialError("worker reported an error: " + (e->is_string() ? e->str : std::string("(non-string)")));
        try {
            const json::Value* rid = r.get("id");
            if (!rid || !rid->is_number() || rid->num != double(id)) throw ExternalPotentialError("reply id does not match request id " + std::to_string(id));
            EvalResult R;
            R.energy = number(r, "energy");
            R.forces = vectors(r, "forces", N);
            const json::Value* vir = r.get("virial");
            const json::Value* str = r.get("stress");
            if (vir) R.virial = mat3(*vir, "virial");
            else if (str) R.virial = mat3(*str, "stress") * s.volume();
            else throw ExternalPotentialError("reply has neither \"virial\" nor \"stress\"");
            return R;
        } catch (const ExternalPotentialError& e) {
            kill_worker();
            throw ExternalPotentialError(std::string(e.what()) + log_tail());
        }
    }

private:
    std::string cmd_;
    std::vector<std::string> names_;
    ExternalOptions opt_;
    std::string worker_name_ = "?";
    double cutoff_ = 0;
    std::vector<double> masses_;
    mutable pid_t pid_ = -1;
    mutable pid_t pgid_ = -1;    // the worker's process group (= its pid at spawn)
    mutable int fd_ = -1;        // our end of the socketpair (worker stdin+stdout)
    mutable int lifeline_ = -1;  // write end of the supervisor's lifeline pipe (closing it ends the worker's group)
    mutable std::string rbuf_;
    mutable long next_id_ = 1;
    mutable std::atomic<bool> busy_{false};
    long log_offset_ = 0;          // where this instance's section of the log starts

    static std::string default_log_path() {
        static std::atomic<int> counter{0};
        const char* t = std::getenv("TMPDIR");
        std::string dir = (t && *t) ? t : "/tmp";
        return dir + "/atomforge-worker-" + std::to_string(::getpid()) + "-" + std::to_string(++counter) + ".log";
    }
    // "exec cmd" lets the worker replace the shell (so signals and the exit status are the worker's own), but only for a plain
    // command: with shell operators or a leading builtin ("cd dir && ...") the shell must stay.
    static std::string shell_command(const std::string& cmd) {
        const bool plain = cmd.find_first_of(";&|<>()`$\n") == std::string::npos;
        const std::size_t b = cmd.find_first_not_of(" \t");
        const std::string first = b == std::string::npos ? "" : cmd.substr(b, cmd.find_first_of(" \t", b) - b);
        static const char* builtins[] = {"cd", "export", "set", "unset", ".", "source", "eval", "ulimit", "umask", "exec", "alias"};
        bool builtin = false;
        for (const char* w : builtins) builtin = builtin || first == w;
        if (!plain || builtin) return cmd;
        // "OMP_NUM_THREADS=1 python3 ..." : exec cannot take an assignment, env can
        const std::size_t eq = first.find('=');
        const bool assign = eq != std::string::npos && eq > 0 && first.find_first_not_of(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_") == eq && !(first[0] >= '0' && first[0] <= '9');
        return (assign ? "exec env " : "exec ") + cmd;
    }

    // The supervisor: /bin/sh running this script, as the leader of a new process group. fd 3 is the read end of a
    // "lifeline" pipe whose only write end atomforge holds; fd 4 is the worker's end of the socket; $1 is the worker
    // command. It starts the worker and a watcher that reads the lifeline. When atomforge exits for any reason (normal
    // exit, Ctrl-C, timeout, SIGHUP, SIGKILL, a crash) the lifeline reads end of file and the watcher kills the whole
    // group at once. When the worker ends, the supervisor ends the same way (same exit status, or the same signal with
    // core dumps off), so atomforge sees the worker's own status. A shell's background job would get /dev/null as stdin,
    // hence the explicit fd 4 redirections. Being an exec'd shell, it shares no memory, descriptors or signal handlers
    // with the host process.
    static const char* supervisor_script() {
        return "w=$1\n"
               "/bin/sh -c \"$w\" 0<&4 1>&4 3<&- 4<&- &\n"
               "pid=$!\n"
               "( exec 0<&3 1>/dev/null 3<&- 4<&-; cat; kill -KILL 0 ) &\n"
               "watch=$!\n"
               "exec 3<&- 4<&- 0</dev/null 1>/dev/null\n"
               "wait \"$pid\"; s=$?\n"
               "kill \"$watch\" 2>/dev/null\n"
               "if [ \"$s\" -gt 128 ] && [ \"$s\" -lt 160 ]; then ulimit -c 0; trap - $((s - 128)) 2>/dev/null; kill -$((s - 128)) $$; fi\n"
               "exit \"$s\"\n";
    }
    static void set_cloexec(int fd) { ::fcntl(fd, F_SETFD, ::fcntl(fd, F_GETFD) | FD_CLOEXEC); }

    void start() {
        // The log is opened here, so an unusable path is reported as such (not as a missing worker command).
        const int logfd = ::open(opt_.log_path.c_str(), O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC, 0644);
        if (logfd < 0) throw ExternalPotentialError("cannot open worker log '" + opt_.log_path + "': " + std::strerror(errno));
        const std::string head = "=== atomforge ExternalPotential (pid " + std::to_string(::getpid()) + "): " + cmd_ + "\n";
        if (::write(logfd, head.data(), head.size()) < 0) { /* the log is best effort */ }
        log_offset_ = long(::lseek(logfd, 0, SEEK_END));
        int sv[2] = {-1, -1}, life[2] = {-1, -1};
        if (::socketpair(AF_UNIX, SOCK_STREAM, 0, sv) != 0 || ::pipe(life) != 0) {
            const int e = errno;
            for (int f : {sv[0], sv[1], life[0], life[1], logfd}) if (f >= 0) ::close(f);
            throw ExternalPotentialError(std::string("cannot create the worker channel: ") + std::strerror(e));
        }
        for (int f : {sv[0], sv[1], life[0], life[1]}) set_cloexec(f);
#if defined(SO_NOSIGPIPE)
        { const int one = 1; ::setsockopt(sv[0], SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof one); }
#endif
        // Everything the child needs is prepared before fork(): it only makes async-signal-safe calls.
        const std::string w = shell_command(cmd_);
        std::vector<char> script(supervisor_script(), supervisor_script() + std::strlen(supervisor_script()) + 1);
        std::vector<char> wc(w.begin(), w.end()); wc.push_back('\0');
        char a0[] = "sh", a1[] = "-c", a3[] = "atomforge-supervisor";
        char* argv[] = {a0, a1, script.data(), a3, wc.data(), nullptr};
        struct rlimit rl{};
        const long maxfd = (::getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY) ? long(std::min<rlim_t>(rl.rlim_cur, 1 << 20)) : 65536;
        const pid_t pid = ::fork();
        if (pid < 0) {
            const int e = errno;
            for (int f : {sv[0], sv[1], life[0], life[1], logfd}) ::close(f);
            throw ExternalPotentialError(std::string("cannot start worker (fork failed): ") + std::strerror(e));
        }
        if (pid == 0) {
            ::setpgid(0, 0);
            // default signal dispositions and an empty mask (exec keeps ignored signals and the mask otherwise: a host
            // that ignores SIGCHLD would break the supervisor's wait)
            for (int sig = 1; sig < NSIG; ++sig) {
                struct sigaction sa{};
                sa.sa_handler = SIG_DFL;
                sigemptyset(&sa.sa_mask);
                ::sigaction(sig, &sa, nullptr);
            }
            sigset_t none; sigemptyset(&none); ::sigprocmask(SIG_SETMASK, &none, nullptr);
            // fds: 0, 1 /dev/null, 2 the log, 3 the lifeline, 4 the socket; nothing else of the host's survives
            const int hi_log = ::fcntl(logfd, F_DUPFD, 10), hi_life = ::fcntl(life[0], F_DUPFD, 10), hi_sock = ::fcntl(sv[1], F_DUPFD, 10);
            const int devnull = ::open("/dev/null", O_RDWR);
            if (hi_log < 0 || hi_life < 0 || hi_sock < 0 || devnull < 0) ::_exit(127);
            ::dup2(devnull, 0); ::dup2(devnull, 1); ::dup2(hi_log, 2); ::dup2(hi_life, 3); ::dup2(hi_sock, 4);
            for (int f = 0; f <= 4; ++f) ::fcntl(f, F_SETFD, 0);   // keep these across exec
#if defined(__linux__) && defined(SYS_close_range)
            if (::syscall(SYS_close_range, 5u, ~0u, 0u) != 0)
#endif
                for (long f = 5; f < maxfd; ++f) ::close(int(f));
            ::execv("/bin/sh", argv);
            ::_exit(127);
        }
        ::setpgid(pid, pid);   // also from the parent, so the group exists before any kill(-pgid)
        ::close(sv[1]);
        ::close(life[0]);
        ::close(logfd);
        pid_ = pgid_ = pid;
        fd_ = sv[0];
        lifeline_ = life[1];
        ::fcntl(fd_, F_SETFL, ::fcntl(fd_, F_GETFL) | O_NONBLOCK);
    }

    void handshake() {
        std::string q = "{\"op\":\"hello\",\"protocol\":" + std::to_string(external_protocol_version) + ",\"elements\":[";
        for (std::size_t i = 0; i < names_.size(); ++i) q += (i ? "," : "") + json::quote(names_[i]);
        q += "]}\n";
        const json::Value r = request(q, opt_.hello_timeout_s, "hello");   // transport errors already carry the log tail
        try {
            if (const json::Value* e = r.get("error"); e && e->type != json::Value::Type::Null) throw ExternalPotentialError("worker refused the handshake: " + (e->is_string() ? e->str : std::string("?")));
            const double pv = number(r, "protocol");
            if (pv != external_protocol_version)
                throw ExternalPotentialError("worker speaks protocol " + json::num(pv) + ", expected " + std::to_string(external_protocol_version));
            if (const json::Value* n = r.get("name"); n && n->is_string()) worker_name_ = n->str;
            if (r.get("cutoff")) cutoff_ = number(r, "cutoff");
            masses_.assign(names_.size(), 0.0);
            for (std::size_t i = 0; i < names_.size(); ++i)
                if (const auto* p = ZhouEAM::lookup(names_[i])) masses_[i] = p->mass;
            if (const json::Value* m = r.get("masses")) {
                if (!m->is_array() || m->arr.size() != names_.size()) throw ExternalPotentialError("\"masses\" must be an array with one entry per element");
                for (std::size_t i = 0; i < names_.size(); ++i) {
                    if (!m->arr[i].is_number() || !std::isfinite(m->arr[i].num)) throw ExternalPotentialError("\"masses\" has a non-finite entry");
                    masses_[i] = m->arr[i].num;
                }
            }
            if (!opt_.masses.empty()) {
                if (opt_.masses.size() != names_.size()) throw ExternalPotentialError("ExternalOptions::masses needs one entry per element");
                masses_ = opt_.masses;
            }
        } catch (const ExternalPotentialError& e) {
            kill_worker();
            throw ExternalPotentialError(std::string(e.what()) + log_tail());
        }
    }

    // Send one line, read one line, parse it. On any transport failure the worker is killed and the error carries the log tail.
    json::Value request(const std::string& line, double timeout_s, const char* what) const {
        using clock = std::chrono::steady_clock;
        // Very large (>= 1e9 s, about 32 years), infinite or NaN timeouts mean "no deadline"; converting them to a
        // clock duration would overflow.
        const bool no_deadline = !(timeout_s < 1e9);
        const auto start = clock::now();
        const auto deadline = no_deadline ? start : start + std::chrono::duration_cast<clock::duration>(std::chrono::duration<double>(timeout_s));
        auto expired = [&]() { return !no_deadline && clock::now() >= deadline; };
        auto remaining_ms = [&]() {   // wake at least once a second to check the child
            if (no_deadline) return 1000;
            const auto left = std::chrono::duration_cast<std::chrono::milliseconds>(deadline - clock::now()).count();
            return left < 0 ? 0 : (left > 1000 ? 1000 : int(left));
        };
        auto die = [&](const std::string& why) {
            const std::string st = reap_status();
            kill_worker();
            throw ExternalPotentialError(std::string(what) + ": " + why + st + log_tail());
        };
        // write
        std::size_t off = 0;
        while (off < line.size()) {
            const ssize_t n = ::send(fd_, line.data() + off, line.size() - off, ATOMFORGE_MSG_NOSIGNAL);
            if (n > 0) { off += std::size_t(n); continue; }
            if (n < 0 && errno == EINTR) continue;
            if (n < 0 && (errno == EAGAIN || errno == EWOULDBLOCK)) {
                if (expired()) die("timed out after " + json::num(timeout_s) + " s writing the request");
                pollfd p{fd_, POLLOUT, 0};
                ::poll(&p, 1, remaining_ms());
                if (worker_exited()) die("worker exited");
                continue;
            }
            if (errno == EPIPE || errno == ECONNRESET) die("worker exited before reading the request");
            die(std::string("cannot write to worker (") + std::strerror(errno) + ")");
        }
        // read one line
        for (;;) {
            const std::size_t nl = rbuf_.find('\n');
            if (nl != std::string::npos) {
                std::string msg = rbuf_.substr(0, nl);
                rbuf_.erase(0, nl + 1);
                if (msg.find_first_not_of(" \t\r") == std::string::npos) continue;   // blank line
                try {
                    json::Value v = json::parse(msg);
                    if (!v.is_object()) throw ExternalPotentialError("reply is not a JSON object");
                    return v;
                } catch (const ExternalPotentialError& e) {
                    die(std::string(e.what()) + ": '" + msg.substr(0, 200) + "'");
                }
            }
            if (rbuf_.size() > (std::size_t(256) << 20)) die("reply line longer than 256 MiB");
            if (expired()) die("timed out after " + json::num(timeout_s) + " s waiting for the reply");
            pollfd p{fd_, POLLIN, 0};
            const int pr = ::poll(&p, 1, remaining_ms());
            if (pr < 0 && errno != EINTR) die(std::string("poll failed (") + std::strerror(errno) + ")");
            if (pr <= 0) {
                // No data. The socket stays open while anything the worker started still holds it, so a worker that exits
                // is noticed here, not by end of file: drain what it wrote, then fail unless that holds a whole reply.
                if (worker_exited()) {
                    for (;;) {
                        char buf[65536];
                        const ssize_t n = ::recv(fd_, buf, sizeof buf, 0);
                        if (n <= 0) break;
                        rbuf_.append(buf, std::size_t(n));
                    }
                    if (rbuf_.find('\n') == std::string::npos) die("worker exited");
                }
                continue;
            }
            char buf[65536];
            const ssize_t n = ::recv(fd_, buf, sizeof buf, 0);
            if (n > 0) { rbuf_.append(buf, std::size_t(n)); continue; }
            if (n == 0 || (n < 0 && errno == ECONNRESET)) die("worker closed its output (exited?)");
            if (errno != EINTR && errno != EAGAIN && errno != EWOULDBLOCK) die(std::string("cannot read from worker (") + std::strerror(errno) + ")");
        }
    }

    // Exit status if the child has exited, as " (exit status N)" / " (killed by signal N)"; "" while it is running.
    mutable std::string exit_note_;
    bool worker_exited() const {
        if (pid_ <= 0) return true;
        int st = 0;
        const pid_t r = ::waitpid(pid_, &st, WNOHANG);
        if (r < 0 && errno == ECHILD) {   // already reaped elsewhere (a host that ignores SIGCHLD): it is gone
            exit_note_ = " (exit status not available)";
            pid_ = -1;
            return true;
        }
        if (r == pid_) {
            if (WIFEXITED(st)) exit_note_ = " (exit status " + std::to_string(WEXITSTATUS(st)) + ")";
            else if (WIFSIGNALED(st)) exit_note_ = " (killed by signal " + std::to_string(WTERMSIG(st)) + ")";
            pid_ = -1;
            return true;
        }
        return false;
    }
    std::string reap_status() const {
        // give a dying worker a moment so the exit status is known
        for (int k = 0; k < 20 && !worker_exited(); ++k) ::usleep(10000);
        return exit_note_;
    }
    void kill_worker() const {
        if (fd_ >= 0) { ::close(fd_); fd_ = -1; }
        if (lifeline_ >= 0) { ::close(lifeline_); lifeline_ = -1; }
        // The whole process group, even if the worker itself has already exited and left children behind.
        if (pgid_ > 0) { ::kill(-pgid_, SIGKILL); pgid_ = -1; }
        if (pid_ > 0) {
            ::kill(pid_, SIGKILL);
            int st = 0;
            ::waitpid(pid_, &st, 0);
            pid_ = -1;
        }
    }
    void shutdown() noexcept {
        if (pid_ > 0 && fd_ >= 0) {
            static const char bye[] = "{\"op\":\"exit\"}\n";
            (void)::send(fd_, bye, sizeof bye - 1, ATOMFORGE_MSG_NOSIGNAL);
            ::shutdown(fd_, SHUT_WR);
            for (int k = 0; k < 200 && !worker_exited(); ++k) ::usleep(10000);   // up to 2 s to exit cleanly
        }
        kill_worker();
    }
    std::string log_tail() const {
        std::ifstream f(opt_.log_path, std::ios::binary);
        if (!f) return "";
        f.seekg(log_offset_);
        std::string all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        if (all.empty()) return "";
        if (all.size() > 1500) all = "..." + all.substr(all.size() - 1500);
        while (!all.empty() && (all.back() == '\n' || all.back() == '\r')) all.pop_back();
        return "\n--- worker stderr (" + opt_.log_path + ") ---\n" + all;
    }

    static double number(const json::Value& r, const char* key) {
        const json::Value* v = r.get(key);
        if (!v) throw ExternalPotentialError(std::string("reply has no \"") + key + "\"");
        if (!v->is_number()) throw ExternalPotentialError(std::string("\"") + key + "\" is not a number");
        if (!std::isfinite(v->num)) throw ExternalPotentialError(std::string("\"") + key + "\" is not finite");
        return v->num;
    }
    static Vec3 vec3(const json::Value& v, const std::string& what) {
        if (!v.is_array() || v.arr.size() != 3) throw ExternalPotentialError(what + " is not a 3-vector");
        Vec3 r;
        for (int k = 0; k < 3; ++k) {
            if (!v.arr[k].is_number()) throw ExternalPotentialError(what + " has a non-number entry");
            if (!std::isfinite(v.arr[k].num)) throw ExternalPotentialError(what + " has a non-finite entry");
            r[k] = v.arr[k].num;
        }
        return r;
    }
    static std::vector<Vec3> vectors(const json::Value& r, const char* key, std::size_t n) {
        const json::Value* v = r.get(key);
        if (!v) throw ExternalPotentialError(std::string("reply has no \"") + key + "\"");
        if (!v->is_array() || v->arr.size() != n)
            throw ExternalPotentialError(std::string("\"") + key + "\" must have " + std::to_string(n) + " entries, got " + (v->is_array() ? std::to_string(v->arr.size()) : std::string("a non-array")));
        std::vector<Vec3> out(n);
        for (std::size_t i = 0; i < n; ++i) out[i] = vec3(v->arr[i], std::string(key) + "[" + std::to_string(i) + "]");
        return out;
    }
    static Mat3 mat3(const json::Value& v, const char* key) {
        if (!v.is_array() || v.arr.size() != 3) throw ExternalPotentialError(std::string("\"") + key + "\" must be a 3x3 array");
        Mat3 m;
        for (int i = 0; i < 3; ++i) {
            const Vec3 row = vec3(v.arr[i], std::string(key) + "[" + std::to_string(i) + "]");
            for (int j = 0; j < 3; ++j) m(i, j) = row[j];
        }
        return m;
    }
};

} // namespace atomforge
