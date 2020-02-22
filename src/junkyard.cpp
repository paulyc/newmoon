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

#include <exception>
#include <cmath>

#define FLOAT128_SUPPORTED 0 // so sayeth clang-darwin-amd64
//#define LONG_DOUBLE_IS_128BIT static_assert(sizeof(long double) == 16))

#if FLOAT128_SUPPORTED
#define float128_t __float128
#else
#define float128_t long double
#endif

/**
 The IEEE 754 standard specifies a binary128 as having:

 Sign bit: 1 bit
 Exponent width: 15 bits
 Significand precision: 113 bits (112 explicitly stored)
 */
struct binary128
{
	binary128() {}
	binary128(int16_t signed_exponent, __uint128_t mantissa)
	{
		__uint128_t mask = 0x0000FFFFFFFFFFFFull;
		mask <<= 64;
		mask |= 0xFFFFFFFFFFFFFFFFull;
		__uint128_t exp_shift = static_cast<__uint128_t>(signed_exponent) << 112;
		irepr = exp_shift | (mantissa & mask);
	}
	union {
		float128_t frepr;
		__uint128_t irepr;
		uint8_t _bytes[16];
		uint16_t signed_exp;
		uint8_t mantissa[14];
	};

	enum ParseState {
		Initial, ReadZero, ReadZeroX, ReadDot, ReadP, Final
	};
	static binary128 from_hexstr(const std::string &s) {
		binary128 val;
		ParseState pstate = Initial;
		std::string buffer;
		for (const char &ch : s) {
			if (pstate == Final) {
				// don't care, ignore
				continue;
			}
			switch (ch) {
				case 'x':
				case 'X':
					if (pstate == ReadZero) {
						pstate = ReadZeroX;
					} else {
						throw std::runtime_error("Parser oops, start with 0x");
					}
					break;
				case '0':
					if (pstate == Initial) {
						pstate = ReadZero;
						continue;
					}
					//yea yea compiler, fall ass through
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				case '7':
				case '8':
				case '9':
				case 'a':
				case 'A':
				case 'b':
				case 'B':
				case 'c':
				case 'C':
				case 'd':
				case 'D':
				case 'e': // nope it's not a special case in a hex string
				case 'E':
				case 'f':
				case 'F':
					if (pstate == ReadZeroX || pstate == ReadDot || pstate == ReadP) {
						buffer += ch;
					}
					break;
				case 'p':
				case 'P':
				case 'L':
				case 'l':
				case '-':
					if (pstate == Initial || pstate == ReadP) {
						buffer += '-'; // put as many as you want idc as long as they're in a row
					} else {
						throw std::runtime_error("Parser oops, misplaced -");
					}
				case '.':
					if (pstate == ReadZeroX) {
						//val.signed_exp = atoi(buffer);
					}
				default:
					throw std::runtime_error("Parser Oops");
			}
		}
		return val;
	}

	binary128& operator=(const char *s) {
		*this = from_hexstr(std::string(s));
		return *this;
	}
};
