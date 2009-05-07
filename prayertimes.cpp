/*-------------------- In the name of God ----------------------*\
 
    PrayerTimes 0.1
    Islamic Prayer Times Calculator
		 
Developed by:
  Mohammad Ebrahim Mohammadi Panah <ebrahim at mohammadi dot ir>

------------------------------------------------------------------

Copyright 2009, Mohammad Ebrahim Mohammadi Panah

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You can get a copy of the GNU General Public License from
http://www.gnu.org/copyleft/gpl.html

\*--------------------------------------------------------------*/

#include "prayertimes.hpp"

#define PROG_NAME "prayertimes"
#define PROG_NAME_FRIENDLY "PrayerTimes"

static const char* TimeName[] =
{
	"Fajr",
	"Sunrise",
	"Dhuhr",
	"Asr",
	"Sunset",
	"Maghrib",
	"Isha",
};

int main()
{
	fprintf(stderr, PROG_NAME_FRIENDLY " %u.%u\n\n", PrayerTimes::VERSION_MAJOR, PrayerTimes::VERSION_MINOR);
	PrayerTimes prayer_times;
	prayer_times.set_calc_method(PrayerTimes::Jafari);
	prayer_times.set_asr_method(PrayerTimes::Shafii);
	prayer_times.set_high_lats_adjust_method(PrayerTimes::MidNight);

	time_t t = time(NULL);
	double times[PrayerTimes::TimesCount];
	double timezone = PrayerTimes::get_effective_timezone(t);
	printf("date     : %s", ctime(&t));
	printf("timezone : %.1lf\n\n", timezone);
	prayer_times.get_prayer_times(t, 35.7061, 51.4358, timezone, times);
	for (int i = 0; i < PrayerTimes::TimesCount; ++i)
		printf("%8s : %s\n", TimeName[i], PrayerTimes::float_time_to_time24(times[i]).c_str());
	return 0;
}
