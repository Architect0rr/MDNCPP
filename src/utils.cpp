#ifndef __MDN_utils__
#define __MDN_utils__

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

#include <filesystem>
// #include <streambuf>
// #include <string>

namespace mdn
{
    enum MPI_TAGS : int
    {
        SANITY,
        STATE
    };

    enum RETURN_CODES : int
    {
        OK,
        ERROR,
        EXIT
    };
}

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
        fs::path outfile = "/";
    };

    struct MC
    {
        spdlog::logger logger = spdlog::logger("default");
        const int rank = 0;
        const int size = 0;
        fs::path cwd = "/";
        const Arguments args;
        MC(const int rank, const int size, const Arguments &args, fs::path cwd) : rank(rank), size(size), args(args), cwd(cwd) {}
        MC() = delete;
        // ~MC(){
        //     delete logger;
        // }
        void setup_logger(const fs::path &log_folder, const int rank)
        {
            fs::path log_file = log_folder / (std::to_string(rank) + ".log");

            auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            console_sink->set_level(spdlog::level::trace);

            auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
            file_sink->set_level(spdlog::level::trace);

            logger = spdlog::logger(std::string("Rank") + std::to_string(rank), {console_sink, file_sink});
            logger.set_level(spdlog::level::trace);
        }
    };
}

#endif // !__MDN_utils__