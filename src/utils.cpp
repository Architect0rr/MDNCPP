#ifndef __MDN_utils__
#define __MDN_utils__

#include "spdlog/spdlog.h"
// #include "spdlog/sinks/stdout_color_sinks.h"
// #include "spdlog/sinks/basic_file_sink.h"

#include <filesystem>
// #include <streambuf>
// #include <string>

namespace mdn::utils
{
    namespace fs = std::filesystem;
    class existing_string_buf : public std::streambuf
    {
    public:
        // Somehow store a pointer to to_append.
        explicit existing_string_buf(std::string &to_append) : m_to_append(&to_append) {}

        virtual int_type overflow(int_type c)
        {
            if (c != EOF)
            {
                m_to_append->push_back(c);
            }
            return c;
        }

        virtual std::streamsize xsputn(const char *s, std::streamsize n)
        {
            m_to_append->insert(m_to_append->end(), s, s + n);
            return n;
        }

    private:
        std::string *m_to_append;
    };

    struct Arguments
    {
        bool debug = false;
        bool trace = false;
        int mode = 0;
        bool verbose = false;
    };

    struct MC
    {
        spdlog::logger *logger = nullptr;
        const int rank = 0;
        const int size = 0;
        fs::path cwd = "/";
        const Arguments args;
        MC(const int rank, const int size, const Arguments &args) : rank(rank), size(size), args(args) {}
        MC() = delete;
        // ~MC(){
        //     delete logger;
        // }
        void add_logger(spdlog::logger *other)
        {
            logger = other;
        }
        void set_cwd(fs::path &pwd)
        {
            cwd = pwd;
        }
    };

    enum MPI_TAGS : int
    {
        SANITY,
        STATE
    };
}

#endif // !__MDN_utils__