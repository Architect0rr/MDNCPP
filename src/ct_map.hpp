#include <iostream>
#include <unordered_map>
#include <vector>
#include <string_view>


enum class CT : int
{
    REQUIRED_BY = -2,
    REQUIRES = -1,
    PID = 0,
    CID = 1,
    MASS = 2,
    VELX = 3,
    VELY = 4,
    VELZ = 5,
    X = 6,
    Y = 7,
    Z = 8,
    PE = 9,
    KE = 10,
    PES = 11,
    KES = 12,
    TEMP = 13,
    TEMPS = 14,
}; // class CT


class CTNodeBase{
public:
    CTNodeBase() = delete;
    constexpr CTNodeBase(const std::string_view& name, std::unordered_map<CT, CTNodeBase*>* mapptr = nullptr): name(name), mapptr(mapptr) {};
    const std::string_view name;

    inline const bool consistent() const {
        return __consistent;
    }

    inline void set(){
        if (!__protected){
            __consistent = true;
        } else {
            #ifdef __MDN__WARNINGS__
                std::cerr << "Attempt to change (set) protected CTNode: " << name << std::endl;
            #endif // __MDN__WARNINGS__
        }
    }

    inline void unset(){
        if (!__protected){
            __consistent = false;
        } else {
            #ifdef __MDN__WARNINGS__
                std::cerr << "Attempt to change (set) protected CTNode: " << name << std::endl;a
            #endif // __MDN__WARNINGS__
        }
    }

    virtual void invalidate() = 0;

    virtual bool validate() = 0;

    inline void set_map(std::unordered_map<CT, CTNodeBase*>* mapptr){
        this->mapptr = mapptr;
    }

    inline void protect(){
        __protected = true;
    }

    inline void unprotect(){
        __protected = false;
    }

protected:
    std::unordered_map<CT, CTNodeBase*>* mapptr = nullptr;
    bool __consistent = false;
    bool __protected = false;

}; // class CTNodeBase


template<CT PROP, CT... Handles>
class CTNode: public CTNodeBase {
public:
    CTNode() = delete;
    constexpr CTNode(const std::string_view& name, std::unordered_map<CT, CTNodeBase*>* mapptr = nullptr): CTNodeBase(name, mapptr) {}
    CTNode(const CTNode& other) = delete;
    CTNode(CTNode&& other) = delete;
    CTNode& operator=(const CTNode& other) = delete;
    CTNode& operator=(CTNode&& other) = delete;
    ~CTNode() = default;

    inline void invalidate() override {
        // __consistent = false;
        invalidate_requires<Handles...>();
    }

    inline bool validate() override {
        validate_requires<Handles...>();
    }

private:

    template<CT ct, CT... Hdls>
    bool validate_requires(){
        if constexpr (ct == CT::REQUIRED_BY){
            return true;
        } else {
            if constexpr (sizeof...(Hdls) == 0){
                if constexpr (ct == CT::REQUIRES) {
                    return true;
                } else {
                    return (*mapptr)[ct]->consistent();
                }
            } else {
                if constexpr (ct == CT::REQUIRES){
                    return validate_requires<Hdls...>();
                } else {
                    return (*mapptr)[ct]->consistent() && validate_requires<Hdls...>();
                }
            }
        }
    }

    template<CT Pr, CT... Hdls>
    void invalidate_requires(){
        if constexpr (sizeof...(Hdls) > 0){
            if constexpr (Pr == CT::REQUIRED_BY){
                invalidate_required<Hdls...>();
            } else {
                if constexpr (Pr == CT::REQUIRES){
                    invalidate_requires<Hdls...>();
                } else {
                    invalidate_requires<Hdls...>();
                }
            }
        }
    }

    template<CT Pr, CT... Hdls>
    requires (Pr != CT::REQUIRED_BY && Pr != CT::REQUIRES)
    void invalidate_required(){
        invalidate_field<Pr>();
        if constexpr (sizeof...(Hdls) > 0)
            invalidate_required<Hdls...>();
    }

    template<CT Pr>
    requires (Pr != CT::REQUIRED_BY && Pr != CT::REQUIRES)
    inline void invalidate_field(){
        if (mapptr == nullptr){
            throw std::runtime_error("Map was not set: mapptr is nullptr");
        }
        (*mapptr)[Pr]->unset();
        (*mapptr)[Pr]->invalidate();
        // mapptr->at(Pr)->invalidate();
    }

}; // class CTNODE


struct CTNodeHandler{
    CTNode<CT::PID,   CT::REQUIRES, CT::REQUIRED_BY>                                         V_PID   = CTNode<CT::PID,   CT::REQUIRES, CT::REQUIRED_BY>                                         ("PID"  );
    CTNode<CT::CID,   CT::REQUIRES, CT::REQUIRED_BY>                                         V_CID   = CTNode<CT::CID,   CT::REQUIRES, CT::REQUIRED_BY>                                         ("CID"  );
    CTNode<CT::MASS,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 V_MASS  = CTNode<CT::MASS,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 ("MASS" );
    CTNode<CT::VELX,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 V_VELX  = CTNode<CT::VELX,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 ("VELX" );
    CTNode<CT::VELY,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 V_VELY  = CTNode<CT::VELY,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 ("VELY" );
    CTNode<CT::VELZ,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 V_VELZ  = CTNode<CT::VELZ,  CT::REQUIRES, CT::REQUIRED_BY, CT::KE>                                 ("VELZ" );
    CTNode<CT::X,     CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 V_X     = CTNode<CT::X,     CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 ("X"    );
    CTNode<CT::Y,     CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 V_Y     = CTNode<CT::Y,     CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 ("Y"    );
    CTNode<CT::Z,     CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 V_Z     = CTNode<CT::Z,     CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 ("Z"    );
    CTNode<CT::PE,    CT::REQUIRES, CT::X, CT::Y, CT::Z, CT::REQUIRED_BY>                    V_PE    = CTNode<CT::PE,    CT::REQUIRES, CT::X, CT::Y, CT::Z, CT::REQUIRED_BY>                    ("PE"   );
    CTNode<CT::KE,    CT::REQUIRES, CT::VELX, CT::VELY, CT::VELZ, CT::REQUIRED_BY, CT::TEMP> V_KE    = CTNode<CT::KE,    CT::REQUIRES, CT::VELX, CT::VELY, CT::VELZ, CT::REQUIRED_BY, CT::TEMP> ("KE"   );
    CTNode<CT::TEMP,  CT::REQUIRES, CT::KE, CT::REQUIRED_BY>                                 V_TEMP  = CTNode<CT::TEMP,  CT::REQUIRES, CT::KE, CT::REQUIRED_BY>                                 ("TEMP" );
    CTNode<CT::TEMPS, CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 V_TEMPS = CTNode<CT::TEMPS, CT::REQUIRES, CT::REQUIRED_BY, CT::PE>                                 ("TEMPS");

    std::unordered_map<CT, CTNodeBase*> map = {
        {CT::PID,   &V_PID},
        {CT::CID,   &V_CID},
        {CT::MASS,  &V_MASS},
        {CT::VELX,  &V_VELX},
        {CT::VELY,  &V_VELY},
        {CT::VELZ,  &V_VELZ},
        {CT::X,     &V_X},
        {CT::Y,     &V_Y},
        {CT::Z,     &V_Z},
        {CT::PE,    &V_PE},
        {CT::KE,    &V_KE},
        {CT::TEMP,  &V_TEMP},
        {CT::TEMPS, &V_TEMPS},
    };

    CTNodeHandler(){
        for (auto [prop, node]: map){
            node->set_map(&map);
        }
    }

};


void print_chain(const std::unordered_map<CT, CTNodeBase*>& nodes){
    for (auto [prop, node]: nodes){
        std::cout << node->name << " " << (node->consistent() ? "true" : "false") << std::endl;
    }
}


int main(int argc, const char** argv) {

    CTNodeHandler handler;

    handler.map[CT::VELX]->set();
    handler.map[CT::VELY]->set();
    handler.map[CT::VELZ]->set();

    handler.map[CT::KE]->set();
    handler.map[CT::TEMP]->set();

    print_chain(handler.map);

    std::cout << "Invalidating VELX" << std::endl;
    handler.map[CT::VELX]->invalidate();
    std::cout << "Invalidated VELX" << std::endl;

    print_chain(handler.map);

    return 0;

} // int main
