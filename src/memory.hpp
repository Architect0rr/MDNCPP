#ifndef __MDN_MEMORY__
#define __MDN_MEMORY__

#include <unordered_set>
#include <unordered_map>
#include <span>

#include "nlohmann/json.hpp"

#include "enums.hpp"

namespace mdn{
    extern spdlog::logger logger;
    using json = nlohmann::json;

    class Memory{
    public:
        static Memory construct(const std::unordered_map<int, PROPS>& _table){
            std::unordered_set<PROPS> avails;
            std::unordered_map<PROPS, int> reverse_table;
            for (const auto& [k, v] : _table) {
                avails.emplace(v);
                reverse_table.emplace(v, k);
            }
            return Memory(_table, reverse_table, avails);
        }

        static std::unordered_map<int, PROPS> parse_props(const json& j) {
            std::unordered_map<int, PROPS> result;

            if (!j.is_object()) throw std::runtime_error("PROPS JSON is not an object");

            for (const auto& [key, value] : j.items())
                result.emplace(std::stoi(key), PROPS(value.get<int>()));
            return result;
        }

        Memory(): valid(false) {}

        Memory(Memory&& o):
            valid(o.valid), _nprops(o._nprops), nlines_allocated(o.nlines_allocated),
            table(std::move(o.table)), reverse_table(std::move(o.reverse_table)), avails(std::move(o.avails)), _data(std::move(o._data)),
            _PID(std::move(o._PID)), _CID(std::move(o._CID)), _MASS(std::move(o._MASS)),
            _VELX(std::move(o._VELX)), _VELY(std::move(o._VELY)), _VELZ(std::move(o._VELZ)),
            _X(std::move(o._X)), _Y(std::move(o._Y)), _Z(std::move(o._Z)),
            _PE(std::move(o._PE)), _KE(std::move(o._KE))
        {
            o.valid = false;
            o._nprops = 0;
            o.nlines_allocated = 0;
        }

        Memory& operator=(Memory&& o){
            valid = o.valid;
            o.valid = false;
            _nprops = o._nprops;
            o._nprops = 0;
            nlines_allocated = o.nlines_allocated;
            o.nlines_allocated = 0;
            table = std::move(o.table);
            reverse_table = std::move(o.reverse_table);
            avails = std::move(o.avails);
            _data = std::move(o._data);
            _PID = std::move(o._PID);
            _CID = std::move(o._CID);
            _MASS = std::move(o._MASS);
            _VELX = std::move(o._VELX);
            _VELY = std::move(o._VELY);
            _VELZ = std::move(o._VELZ);
            _X = std::move(o._X);
            _Y = std::move(o._Y);
            _Z = std::move(o._Z);
            _PE = std::move(o._PE);
            _KE = std::move(o._KE);
            return *this;
        }

        inline double* get() {return _data.get();}
        inline const double* const get() const {return _data.get();}

        inline double* data() {
            if (valid && nlines_allocated) return _data.get();
            else throw std::runtime_error("Invalid object state");
        }
        inline const double* const data() const {
            if (valid && nlines_allocated) return _data.get();
            else throw std::runtime_error("Invalid object state");
        }

        void allocate(const uint64_t& count){
            _data = std::unique_ptr<double[]>(new double[count*_nprops]);
            nlines_allocated = count;
            // for (const auto& [k, v] : reverse_table){

            // }
            for (const auto& [k, v] : table){
                switch(v) {
                    case PROPS::PID:
                        _PID = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::CID:
                        _CID = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::MASS:
                        _MASS = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::VELX:
                        _VELX = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::VELY:
                        _VELY = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::VELZ:
                        _VELZ = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::X:
                        _X = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::Y:
                        _Y = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::Z:
                        _Z = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::PE:
                        _PE = std::span<const double>(_data.get() + k*count, count);
                        break;
                    case PROPS::KE:
                        _KE = std::span<const double>(_data.get() + k*count, count);
                        break;
                    default:
                        throw std::runtime_error("Specified property cannot be parsed");  // usually unreachable code
                }
            }
        }

        inline const std::span<const double>& PID() const {return _PID;}
        inline bool check_PID() {return avails.contains(PROPS::PID);}

        const std::span<const double>& get_PID() const {
            if (avails.contains(PROPS::PID)) return _PID;
            else throw std::runtime_error("There is no PID property");
        }

        inline const std::span<const double>& CID() const {return _CID;}
        inline bool check_CID() {return avails.contains(PROPS::CID);}

        const std::span<const double>& get_CID() const {
            if (avails.contains(PROPS::CID)) return _CID;
            else throw std::runtime_error("There is no CID property");
        }

        inline const std::span<const double>& MASS() const {return _MASS;}
        inline bool check_MASS() {return avails.contains(PROPS::MASS);}

        const std::span<const double>& get_MASS() const {
            if (avails.contains(PROPS::MASS)) return _MASS;
            else throw std::runtime_error("There is no MASS property");
        }

        inline const std::span<const double>& VELX() const {return _VELX;}
        inline bool check_VELX() {return avails.contains(PROPS::VELX);}

        const std::span<const double>& get_VELX() const {
            if (avails.contains(PROPS::VELX)) return _VELX;
            else throw std::runtime_error("There is no VELX property");
        }

        inline const std::span<const double>& VELY() const {return _VELY;}
        inline bool check_VELY() {return avails.contains(PROPS::VELY);}

        const std::span<const double>& get_VELY() const {
            if (avails.contains(PROPS::VELY)) return _VELY;
            else throw std::runtime_error("There is no VELY property");
        }

        inline const std::span<const double>& VELZ() const {return _VELZ;}
        inline bool check_VELZ() {return avails.contains(PROPS::VELZ);}

        const std::span<const double>& get_VELZ() const {
            if (avails.contains(PROPS::VELZ)) return _VELZ;
            else throw std::runtime_error("There is no VELZ property");
        }

        inline const std::span<const double>& X() const  {return _X;}
        inline bool check_X() {return avails.contains(PROPS::X);}

        const std::span<const double>& get_X() const {
            if (avails.contains(PROPS::X)) return _X;
            else throw std::runtime_error("There is no X property");
        }

        inline const std::span<const double>& Y() const {return _Y;}
        inline bool check_Y() {return avails.contains(PROPS::Y);}

        const std::span<const double>& get_Y() const {
            if (avails.contains(PROPS::Y)) return _Y;
            else throw std::runtime_error("There is no Y property");
        }

        inline const std::span<const double>& Z() const {return _Z;}
        inline bool check_Z() {return avails.contains(PROPS::Z);}

        const std::span<const double>& get_Z() const {
            if (avails.contains(PROPS::Z)) return _Z;
            else throw std::runtime_error("There is no Z property");
        }

        inline const std::span<const double>& PE() const {return _PE;}
        inline bool check_PE() {return avails.contains(PROPS::PE);}

        const std::span<const double>& get_PE() const {
            if (avails.contains(PROPS::PE)) return _PE;
            else throw std::runtime_error("There is no PE property");
        }

        inline const std::span<const double>& KE() const {return _KE;}
        inline bool check_KE() {return avails.contains(PROPS::KE);}

        const std::span<const double>& get_KE() const {
            if (avails.contains(PROPS::KE)) return _KE;
            else throw std::runtime_error("There is no KE property");
        }

        inline const int nprops() {return _nprops;}

    private:
        Memory(const std::unordered_map<int, PROPS>& _table, const std::unordered_map<PROPS, int>& _reverse_table, const std::unordered_set<PROPS>& _avails):
            table(_table), reverse_table(_reverse_table), _nprops(_table.size()), avails(_avails), valid(true) { }
        Memory(const Memory&) = delete;
        Memory& operator=(const Memory&) = delete;

        bool valid = false;
        int _nprops = 0, nlines_allocated = 0;
        std::unordered_map<int, PROPS> table;
        std::unordered_map<PROPS, int> reverse_table;
        std::unordered_set<PROPS> avails;

        std::unique_ptr<double[]> _data;

        std::span<const double> _PID;
        std::span<const double> _CID;
        std::span<const double> _MASS;
        std::span<const double> _VELX;
        std::span<const double> _VELY;
        std::span<const double> _VELZ;
        std::span<const double> _X;
        std::span<const double> _Y;
        std::span<const double> _Z;
        std::span<const double> _PE;
        std::span<const double> _KE;
    };

}  // namepace mdn

#endif // !__MDN_MEMORY__
