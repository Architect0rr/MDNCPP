#ifndef __MDN_enums_HPP__
#define __MDN_enums_HPP__

namespace mdn{
    enum MPI_TAGS : int
    {
        SANITY,
        STATE,
        STOR
    };

    enum RETURN_CODES : int
    {
        OK,
        ERROR,
        EXIT
    };
}

#endif // !__MDN_enums_HPP__