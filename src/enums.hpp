#ifndef __MDN_enums_HPP__
#define __MDN_enums_HPP__

namespace mdn{
    enum class MPI_TAGS : int
    {
        SANITY = 0,
        STATE,
        STOR,
        OUTFILE_EXCHANGE
    };

    enum class RETURN_CODES : int
    {
        OK = 0,
        ERROR = 1,
        EXIT = 2
    };

    enum class PROPS : int
    {
        REQUIRED_BY = -2,
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
    };

}

#endif // !__MDN_enums_HPP__