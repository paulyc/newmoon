/**
 * part of newmoon, moon phase calculator
 *
 * Copyright (C) 2020 Paul Ciarlo <paul.ciarlo@gmail.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 **/

#ifndef PAULYC_JD_CLOCK_HPP
#define PAULYC_JD_CLOCK_HPP

#include <chrono>
#include <vector>
#include <ratio>
#include <sstream>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <cstring>

struct jd_clock : public std::chrono::steady_clock
{
    typedef std::chrono::duration<long double, std::ratio<1,86400>> duration;
    typedef duration::rep                                 rep;
    typedef duration::period                              period;
    typedef std::chrono::time_point<jd_clock, duration> 	time_point;

    /*
     * delta-T =
          ET - UT   prior to 1984
          TDT - UT  1984 - 2000
          TT - UT   from 2001 and on

          delta-UT =  UT - UTC
          unixtime: no leap seconds, so UT
          */
    // epoch = noon Universal Time (UT) Monday, 1 January 4713 BCE (-4712) Julian
    // J = 2000 + (Julian date − 2451545.0) ÷ 365.25
    // jd2000 epoch = 2451545.0 = TT 12:00:00.000 (UTC 11:58:55.816) 1 January 2000 CE Gregorian
    /*
        The Gregorian date January 1, 2000, at 12:00 TT (Terrestrial Time).
        The Julian date 2451545.0 TT (Terrestrial Time).[12]
        January 1, 2000, 11:59:27.816 TAI (International Atomic Time).[13]
        January 1, 2000, 11:58:55.816 UTC (Coordinated Universal Time).[14]
        Smoothed historical measurements of ΔT using total solar eclipses are about
        +17190 s in the year −500 (501 BC),
        +10580 s in 0 (1 BC),
        +5710 s in 500,
        +1570 s in 1000, and
        +200  s in 1500.
        After the invention of the telescope, measurements were made by observing occultations of stars by the Moon,
        which allowed the derivation of more closely spaced and more accurate values for ΔT.
        ΔT continued to decrease until it reached a plateau of +11 ± 6 s between 1680 and 1866.
        For about three decades immediately before 1902 it was negative, reaching -6.64 s. Then it increased to
        +63.83 s in January 2000
        +68.97 s in January 2018.
    */
    static constexpr bool is_steady = true;
    static constexpr long double SECONDS_PER_JDAY = 86400.0l;
    static constexpr long double JDAYS_PER_YEAR = 365.25l;
    static constexpr long double SECONDS_PER_JYEAR = SECONDS_PER_JDAY * JDAYS_PER_YEAR;
    static constexpr long double JD2000_EPOCH_JD = 2451545.0l;
    static constexpr long JD2000_EPOCH_JDS = 211813488000l;
    static constexpr long double JD2000_EPOCH_UNIXTIME = 946727935.816l;//946728000l;//
    static constexpr long double UNIX_EPOCH_JD = JD2000_EPOCH_JD - JD2000_EPOCH_UNIXTIME/SECONDS_PER_JDAY;
    static constexpr long double JD2000_DELTA_T = 63.83l;
    static constexpr long double DELTA_T_2018 = 68.97l;
    static constexpr long double leap_seconds_since_unix_epoch = 27.0l;

    static time_point now() noexcept {
        return from_unixtime(std::chrono::system_clock::now());
    }

    static time_point from_unixtime(const std::chrono::system_clock::time_point &t) {
        std::time_t unixtime = std::chrono::system_clock::to_time_t(t);
        return time_point(duration(UNIX_EPOCH_JD + unixtime/SECONDS_PER_JDAY));
    }

    struct YearDT {
        double year;
        double dt;
    };

    static const std::vector<YearDT> DELTA_T_YR;

    static double delta_t_lerp(double year) {
        size_t point = 1;
        while (point < jd_clock::DELTA_T_YR.size() - 1 && year > jd_clock::DELTA_T_YR[point].year) {
            ++point;
        }
        //std::cerr << "year " << year << std::endl;
        const double slope = (jd_clock::DELTA_T_YR[point].dt - jd_clock::DELTA_T_YR[point-1].dt) / (jd_clock::DELTA_T_YR[point].year - jd_clock::DELTA_T_YR[point-1].year);
        //std::cerr << "slope " << slope << std::endl;
        const double dy = jd_clock::DELTA_T_YR[point].year - year;
        //std::cerr << "dy " << dy << std::endl;
        const double dt = jd_clock::DELTA_T_YR[point].dt - dy * slope;
        //std::cerr << "dt " << dt << std::endl;
        return dt;
    }
    static void test_delta_t_lerp() {
        assert(delta_t_lerp(-700) > 17190);
        assert(delta_t_lerp(1250) > 200 && delta_t_lerp(1250) < 1570);
        assert(delta_t_lerp(2030) > 68.97);
    }
};

inline static std::ostream& operator<<(std::ostream &os, const std::chrono::system_clock::time_point &rhs)
{
    const std::time_t tt = std::chrono::system_clock::to_time_t(rhs);
    os << std::put_time(std::gmtime(&tt),"%FT%T%z (%Z)");
    return os;
}

inline static std::ostream& operator<<(std::ostream &os, const jd_clock::time_point &rhs)
{
    os << std::fixed << std::setw( 11 ) << std::setprecision( 6 )
       << std::setfill( '0' ) << rhs.time_since_epoch().count();
    return os;
}

inline static std::string string_format(const char *format, ...)
{
    char buf[1024];
    va_list args;
    va_start(args, format);
    int len = std::vsnprintf(buf, sizeof(buf), format, args);
    va_end(args);

    if (len < 0) {
        throw std::runtime_error(strerror(errno));
    } else if (static_cast<size_t>(len) >= sizeof(buf)) {
        std::string formatted(static_cast<size_t>(len) + 1, '\0');
        va_start(args, format);
        std::vsnprintf(const_cast<char*>(formatted.data()), static_cast<size_t>(len) + 1, format, args);
        va_end(args);
        return formatted;
    } else {
        return std::string(buf);
    }
}

#endif /* PAULYC_JD_CLOCK_HPP */
