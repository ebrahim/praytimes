/*---------------------------------------------------------*\
 
      PrayerTimes: Islamic Prayer Times Calculator

-------------------------------------------------------------

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

\*---------------------------------------------------------*/

#include "prayertimes.hpp"

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
	PrayerTimes prayer_times;
	prayer_times.set_calc_method(PrayerTimes::Jafari);
	prayer_times.set_asr_method(PrayerTimes::Shafii);
	prayer_times.set_high_lats_adjust_method(PrayerTimes::MidNight);

	time_t t = time(NULL);
	double times[PrayerTimes::TimesCount];
	prayer_times.get_prayer_times(t, 35.7061, 51.4358, 4.5, times);
	for (int i = 0; i < PrayerTimes::TimesCount; ++i)
		printf("%8s : %10lf : %s\n", TimeName[i], times[i], PrayerTimes::float_to_time24(times[i]).c_str());
	return 0;
}
