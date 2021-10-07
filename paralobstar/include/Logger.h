//
// Created by Johannes Martin on 27.09.21.
//

#ifndef PARALOBSTAR_LOGGER_H
#define PARALOBSTAR_LOGGER_H

#include <iostream>
#include <string>

#include "global.h"

namespace Color {

    enum Code {
        FG_DEFAULT = 39,
        FG_BLACK = 30,
        FG_RED = 31,
        FG_GREEN = 32,
        FG_YELLOW = 33,
        FG_BLUE = 34,
        FG_MAGENTA = 35,
        FG_CYAN = 36,
        FG_LIGHT_GRAY = 37,
        FG_DARK_GRAY = 90,
        FG_LIGHT_RED = 91,
        FG_LIGHT_GREEN = 92,
        FG_LIGHT_YELLOW = 93,
        FG_LIGHT_BLUE = 94,
        FG_LIGHT_MAGENTA = 95,
        FG_LIGHT_CYAN = 96,
        FG_WHITE = 97,
        BG_RED = 41,
        BG_GREEN = 42,
        BG_BLUE = 44,
        BG_DEFAULT = 49
    };

    class Modifier {
    public:
        Code code;
        Modifier(Code pCode);
        friend std::ostream& operator<<(std::ostream& os, const Color::Modifier& mod);
    };
}

enum typelog {
    DEBUG,
    INFO,
    WARN,
    ERROR
};

struct structlog {
    bool headers = false;
    typelog level = WARN;
    int myRank = -1; // don't use MPI by default
    int outputRank = -1;
};

extern structlog LOGCFG;

class Logger {
public:
    Logger() {}
    Logger(typelog type);
    ~Logger();

    template<class T> Logger &operator<<(const T &msg) {
        if (msglevel >= LOGCFG.level && (LOGCFG.myRank == LOGCFG.outputRank || LOGCFG.outputRank == -1)) {
            std::cout << msg;
            opened = true;
        }
        return *this;
    }

    Logger &operator<<(const std::uint_fast64_t &key) {
        int level = global::maxTreeLvl;
        if (msglevel >= LOGCFG.level && (LOGCFG.myRank == LOGCFG.outputRank || LOGCFG.outputRank == -1)) {
            int levels [level];
            for (int i = 0; i<level; i++) {
                levels[i] = (key >> 3*i) & (unsigned long)7;
            }
            std::string msg = "#|";
            for (int i = level-1; i>=0; i--) {
                msg += std::to_string(levels[i]);
                msg += "|";
            }
            std::cout << msg;
            opened = true;
        }
        return *this;
    }



private:
    bool opened = false;
    typelog msglevel = DEBUG;
    inline std::string getLabel(typelog type);
    inline Color::Modifier getColor(typelog type);
};


#endif //PARALOBSTAR_LOGGER_H
