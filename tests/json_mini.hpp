// Minimal strict JSON parser for tests only (objects keep insertion order). Throws std::runtime_error on bad input.
#pragma once
#include <cstdlib>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mini_json {

struct Value {
    enum Type { Null, Bool, Number, String, Array, Object } type = Null;
    bool b = false;
    double n = 0;
    std::string s;
    std::vector<Value> a;
    std::vector<std::pair<std::string, Value>> o;

    bool is_null() const { return type == Null; }
    bool has(const std::string& k) const { for (auto& kv : o) if (kv.first == k) return true; return false; }
    const Value& at(const std::string& k) const {
        for (auto& kv : o) if (kv.first == k) return kv.second;
        throw std::runtime_error("missing key: " + k);
    }
    // dotted-path lookup through objects, e.g. "quantities.Q_dielectric.value"
    const Value& path(const std::string& p) const {
        const Value* v = this; std::size_t i = 0;
        while (i <= p.size()) {
            const std::size_t j = p.find('.', i);
            const std::string k = p.substr(i, j == std::string::npos ? std::string::npos : j - i);
            v = &v->at(k);
            if (j == std::string::npos) break;
            i = j + 1;
        }
        return *v;
    }
};

class Parser {
public:
    explicit Parser(const std::string& t) : t_(t) {}
    Value parse() { Value v = value(); ws(); if (p_ != t_.size()) fail("trailing characters"); return v; }
private:
    const std::string& t_; std::size_t p_ = 0;
    [[noreturn]] void fail(const std::string& m) const { throw std::runtime_error("json: " + m + " at offset " + std::to_string(p_)); }
    void ws() { while (p_ < t_.size() && (t_[p_] == ' ' || t_[p_] == '\n' || t_[p_] == '\t' || t_[p_] == '\r')) ++p_; }
    char peek() { ws(); if (p_ >= t_.size()) fail("unexpected end"); return t_[p_]; }
    void expect(char c) { if (peek() != c) fail(std::string("expected '") + c + "'"); ++p_; }
    bool lit(const char* w) { std::string s(w); if (t_.compare(p_, s.size(), s) == 0) { p_ += s.size(); return true; } return false; }
    Value value() {
        const char c = peek(); Value v;
        if (c == '{') {
            ++p_; v.type = Value::Object;
            if (peek() == '}') { ++p_; return v; }
            for (;;) {
                if (peek() != '"') fail("object key must be a string");
                std::string k = str(); expect(':');
                v.o.emplace_back(std::move(k), value());
                if (peek() == ',') { ++p_; continue; }
                expect('}'); return v;
            }
        }
        if (c == '[') {
            ++p_; v.type = Value::Array;
            if (peek() == ']') { ++p_; return v; }
            for (;;) { v.a.push_back(value()); if (peek() == ',') { ++p_; continue; } expect(']'); return v; }
        }
        if (c == '"') { v.type = Value::String; v.s = str(); return v; }
        if (lit("true")) { v.type = Value::Bool; v.b = true; return v; }
        if (lit("false")) { v.type = Value::Bool; v.b = false; return v; }
        if (lit("null")) { return v; }
        const char* b = t_.c_str() + p_; char* e = nullptr;
        if (!((c >= '0' && c <= '9') || c == '-')) fail("unexpected character");
        v.n = std::strtod(b, &e);
        if (e == b) fail("bad number");
        p_ += std::size_t(e - b); v.type = Value::Number; return v;
    }
    std::string str() {
        expect('"'); std::string r;
        while (p_ < t_.size() && t_[p_] != '"') {
            char c = t_[p_++];
            if (c == '\\') {
                if (p_ >= t_.size()) fail("bad escape");
                char e = t_[p_++];
                switch (e) {
                    case '"': r += '"'; break; case '\\': r += '\\'; break; case '/': r += '/'; break;
                    case 'n': r += '\n'; break; case 't': r += '\t'; break; case 'r': r += '\r'; break;
                    case 'b': r += '\b'; break; case 'f': r += '\f'; break;
                    case 'u': {
                        if (p_ + 4 > t_.size()) fail("bad \\u escape");
                        const unsigned cp = unsigned(std::strtoul(t_.substr(p_, 4).c_str(), nullptr, 16)); p_ += 4;
                        if (cp < 0x80) r += char(cp);
                        else if (cp < 0x800) { r += char(0xC0 | (cp >> 6)); r += char(0x80 | (cp & 0x3F)); }
                        else { r += char(0xE0 | (cp >> 12)); r += char(0x80 | ((cp >> 6) & 0x3F)); r += char(0x80 | (cp & 0x3F)); }
                        break;
                    }
                    default: fail("bad escape");
                }
            } else if (static_cast<unsigned char>(c) < 0x20) fail("raw control character in string");
            else r += c;
        }
        if (p_ >= t_.size()) fail("unterminated string");
        ++p_; return r;
    }
};

inline Value parse(const std::string& text) { return Parser(text).parse(); }

} // namespace mini_json
