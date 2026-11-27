// atomforge — ExternalPotential: energies, forces and virials from a persistent worker process
//
// The worker is spawned once (through /bin/sh -c "exec <command>") and talks over its stdin/stdout, one JSON object per line.
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
// did not converge should not end a search). POSIX only (posix_spawn, socketpair, poll). Not thread-safe.
#pragma once
#include "potential.hpp"
#include <cerrno>
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
#include <spawn.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

extern char** environ;

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
            const std::size_t start = i_;
            while (i_ < s_.size() && std::strchr("+-.eE0123456789", s_[i_])) ++i_;
            const std::string tok = s_.substr(start, i_ - start);
            char* end = nullptr;
            v.type = Value::Type::Number;
            v.num = std::strtod(tok.c_str(), &end);
            if (!end || *end != '\0') fail("bad number '" + tok + "'");
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
inline std::string num(double x) { char b[32]; std::snprintf(b, sizeof b, "%.17g", x); return b; }

} // namespace json

// ---------------------------------------------------------------- the potential
struct ExternalOptions {
    double hello_timeout_s = 60;     // start-up and handshake (imports of PySCF / torch can be slow)
    double eval_timeout_s = 3600;    // one evaluate request; DFT evaluations can take minutes
    std::string log_path;            // worker stderr; empty -> <tmp>/atomforge-worker-<pid>.log
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

    EvalResult evaluate(const Structure& s, NeighborList&) const override {
        if (pid_ <= 0) throw ExternalPotentialError("external worker is not running (it failed earlier; see " + opt_.log_path + ")");
        const std::size_t N = s.size();
        if (s.species.size() != N) throw ExternalPotentialError("structure has " + std::to_string(s.species.size()) + " species for " + std::to_string(N) + " atoms");
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
        if (const json::Value* e = r.get("error"); e && r.get("id") && r.get("id")->is_number() && r.get("id")->num == double(id))
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
    mutable int fd_ = -1;        // our end of the socketpair (worker stdin+stdout)
    mutable std::string rbuf_;
    mutable long next_id_ = 1;

    static std::string default_log_path() {
        const char* t = std::getenv("TMPDIR");
        std::string dir = (t && *t) ? t : "/tmp";
        return dir + "/atomforge-worker-" + std::to_string(::getpid()) + ".log";
    }

    void start() {
        int sv[2];
        if (::socketpair(AF_UNIX, SOCK_STREAM | SOCK_CLOEXEC, 0, sv) != 0)
            throw ExternalPotentialError(std::string("socketpair failed: ") + std::strerror(errno));
        posix_spawn_file_actions_t fa;
        posix_spawn_file_actions_init(&fa);
        posix_spawn_file_actions_adddup2(&fa, sv[1], STDIN_FILENO);
        posix_spawn_file_actions_adddup2(&fa, sv[1], STDOUT_FILENO);
        posix_spawn_file_actions_addopen(&fa, STDERR_FILENO, opt_.log_path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);
        const std::string shcmd = "exec " + cmd_;
        char a0[] = "sh", a1[] = "-c";
        std::vector<char> c(shcmd.begin(), shcmd.end()); c.push_back('\0');
        char* argv[] = {a0, a1, c.data(), nullptr};
        pid_t pid = -1;
        const int rc = ::posix_spawn(&pid, "/bin/sh", &fa, nullptr, argv, environ);
        posix_spawn_file_actions_destroy(&fa);
        ::close(sv[1]);
        if (rc != 0) {
            ::close(sv[0]);
            throw ExternalPotentialError("cannot start worker '" + cmd_ + "': " + std::strerror(rc));
        }
        pid_ = pid;
        fd_ = sv[0];
        ::fcntl(fd_, F_SETFL, ::fcntl(fd_, F_GETFL) | O_NONBLOCK);
    }

    void handshake() {
        std::string q = "{\"op\":\"hello\",\"protocol\":" + std::to_string(external_protocol_version) + ",\"elements\":[";
        for (std::size_t i = 0; i < names_.size(); ++i) q += (i ? "," : "") + json::quote(names_[i]);
        q += "]}\n";
        const json::Value r = request(q, opt_.hello_timeout_s, "hello");   // transport errors already carry the log tail
        try {
            if (const json::Value* e = r.get("error")) throw ExternalPotentialError("worker refused the handshake: " + (e->is_string() ? e->str : std::string("?")));
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
        const auto deadline = clock::now() + std::chrono::duration_cast<clock::duration>(std::chrono::duration<double>(timeout_s));
        auto remaining_ms = [&]() {
            const auto left = std::chrono::duration_cast<std::chrono::milliseconds>(deadline - clock::now()).count();
            return left < 0 ? 0 : (left > 1000 ? 1000 : int(left));   // wake at least once a second to check the child
        };
        auto die = [&](const std::string& why) {
            const std::string st = reap_status();
            kill_worker();
            throw ExternalPotentialError(std::string(what) + ": " + why + st + log_tail());
        };
        // write
        std::size_t off = 0;
        while (off < line.size()) {
            const ssize_t n = ::send(fd_, line.data() + off, line.size() - off, MSG_NOSIGNAL);
            if (n > 0) { off += std::size_t(n); continue; }
            if (n < 0 && errno == EINTR) continue;
            if (n < 0 && (errno == EAGAIN || errno == EWOULDBLOCK)) {
                if (clock::now() >= deadline) die("timed out after " + json::num(timeout_s) + " s writing the request");
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
            if (clock::now() >= deadline) die("timed out after " + json::num(timeout_s) + " s waiting for the reply");
            pollfd p{fd_, POLLIN, 0};
            const int pr = ::poll(&p, 1, remaining_ms());
            if (pr < 0 && errno != EINTR) die(std::string("poll failed (") + std::strerror(errno) + ")");
            if (pr <= 0) continue;
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
            (void)::send(fd_, bye, sizeof bye - 1, MSG_NOSIGNAL);
            ::shutdown(fd_, SHUT_WR);
            for (int k = 0; k < 200 && !worker_exited(); ++k) ::usleep(10000);   // up to 2 s to exit cleanly
        }
        kill_worker();
    }
    std::string log_tail() const {
        std::ifstream f(opt_.log_path, std::ios::binary);
        if (!f) return "";
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
