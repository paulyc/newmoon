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
#include <memory>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <cstring>

struct jd_clock
{
    typedef long double                         rep;
    typedef std::ratio<86400>                   period;
    typedef std::chrono::duration<rep, period>  duration;
    typedef std::chrono::time_point<jd_clock>   time_point;

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
/*
    static constexpr std::time_t jd2000_unixtime() {
        std::tm tm_jd2000 = {
                .tm_sec = 56,
                .tm_min = 58,
                .tm_hour = 23,
                .tm_mday = 1,
                .tm_mon = 0,
                .tm_year = 2000,
                .tm_isdst = 0,
                .tm_gmtoff = 0,
                .tm_zone = "UTC",
            };
        return mktime(&tm_jd2000);
    }
*/
    static constexpr long double SECONDS_PER_JYEAR = SECONDS_PER_JDAY * JDAYS_PER_YEAR;
    static constexpr long double JD2000_EPOCH_JD = 2451545.0l;
    static constexpr long JD2000_EPOCH_JDS = 211813488000l;
    static constexpr long double JD2000_EPOCH_UNIXTIME = 946727935.816l;//946728000l;//
    static constexpr long double UNIX_EPOCH_JD = JD2000_EPOCH_JD - JD2000_EPOCH_UNIXTIME/SECONDS_PER_JDAY;
    static constexpr long double JD2000_DELTA_T = 63.83l;
    static constexpr long double DELTA_T_2018 = 68.97l;
    static constexpr long double leap_seconds_since_unix_epoch = 27.0l;

    static time_point now() noexcept {
        return from_system_clock(std::chrono::system_clock::now());
    }

    static time_point from_system_clock(const std::chrono::system_clock::time_point &t) {
        return from_time_t(std::chrono::system_clock::to_time_t(t));
    }

    static time_point from_time_t(std::time_t t) {
        return time_point(duration(UNIX_EPOCH_JD + t/jd_clock::SECONDS_PER_JDAY));
    }

    static std::time_t to_time_t(jd_clock::time_point &jd) {
        return static_cast<std::time_t>((jd.time_since_epoch().count() - UNIX_EPOCH_JD) * SECONDS_PER_JDAY);
    }

    static std::chrono::system_clock::time_point to_system_clock(jd_clock::time_point &jd) {
        return std::chrono::system_clock::from_time_t(to_time_t(jd));
    }

    struct YearDT {
        double year;
        double dt;
    };

    static constexpr jd_clock::YearDT DELTA_T_YR[] = {
        {-500, 17190},
        {0, 10580},
        {500,5710},
        {1000,1570},
        {1500,200},
        {1900,-6.64},
        {2000,63.83},
        {2018,68.97}
    };

    // probably broken
    static double delta_t_lerp(double year) {
        size_t point = 1;
        while (point < sizeof(DELTA_T_YR) - 1 && year > DELTA_T_YR[point].year) {
            ++point;
        }
        //std::cerr << "year " << year << std::endl;
        const double slope = (DELTA_T_YR[point].dt - DELTA_T_YR[point-1].dt) / (DELTA_T_YR[point].year - DELTA_T_YR[point-1].year);
        //std::cerr << "slope " << slope << std::endl;
        const double dy = DELTA_T_YR[point].year - year;
        //std::cerr << "dy " << dy << std::endl;
        const double dt = DELTA_T_YR[point].dt - dy * slope;
        //std::cerr << "dt " << dt << std::endl;
        return dt;
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

struct string_formatter
{
    std::string _fmt;

    explicit string_formatter(std::string fmt) : _fmt(fmt) {}

    template <typename ... TArgs>
    std::string format(TArgs ... args) {
        std::size_t buflen = _fmt.length() + 1024;
        std::unique_ptr<char[]> buf = std::make_unique<char[]>(buflen);

       // int len = std::snprintf(buf.get(), buflen, "%s", _fmt.c_str());
        int len = std::snprintf(buf.get(), buflen, _fmt.c_str(), args...);

        if (len < 0) {
            throw std::runtime_error(strerror(errno));
        } else if (static_cast<size_t>(len) >= buflen) {
            buflen = static_cast<std::size_t>(len) + 1;
            buf = std::make_unique<char[]>(buflen);
            std::snprintf(buf.get(), buflen, _fmt.c_str(), args...);
        }
        return std::string(buf.get());
    }
};

template <typename T, T ... Fmt>
constexpr auto operator""_fmt()
{
    return string_formatter(std::string({Fmt...}));
}

#endif /* PAULYC_JD_CLOCK_HPP */
