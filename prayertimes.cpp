/*-------------------- In the name of God ----------------------*\

    PrayerTimes 1.1
    Islamic prayer times calculator

Developed by:
  Ebrahim Mohammadi <ebrahim at mohammadi dot ir>

------------------------------------------------------------------

Copyright 2011, Ebrahim Mohammadi

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

#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <getopt.h>

#include "prayertimes.hpp"

#define PROG_NAME "prayertimes"
#define PROG_NAME_FRIENDLY "PrayerTimes"
#define PROG_VERSION "1.1"

using prayertimes::PrayerTimes;

static const char* const TimeName[] =
{
	"Imsak",
	"Fajr",
	"Sunrise",
	"Dhuhr",
	"Asr",
	"Sunset",
	"Maghrib",
	"Isha",
	"Midnight",
};

static const char* const CalculationMethodName[] =
{
	"MWL",
	"ISNA",
	"Egypt",
	"Makkah",
	"Karachi",
	"Jafari",
	"Tehran",
	"Custom",
};

void print_help(FILE* f)
{
	fputs(PROG_NAME_FRIENDLY " " PROG_VERSION "\n\n", stderr);
	fputs("Usage: " PROG_NAME " options...\n"
	      "\n"
	      " Options\n"
	      "    --help                      -h  you're reading it\n"
	      "    --version                   -v  prints name and version, then exits\n"
	      "    --date arg                  -d  get prayer times for arbitrary date\n"
	      "    --timezone arg              -z  get prayer times for arbitrary timezone\n"
	      "  * --latitude arg              -l  latitude of desired location\n"
	      "  * --longitude arg             -n  longitude of desired location\n"
	      "    --elevation arg             -e  elevation of desired location\n"
	      "    --calc-method arg           -c  select prayer time calculation method\n"
	      "    --asr-juristics-method arg  -a  select Juristic method for calculating Asr prayer time\n"
	      "    --high-lats-method arg      -i  select adjusting method for higher latitude\n"
	      " ** --imsak-minutes arg             minutes before Fajr for calculating Imsak time\n"
	      "    --dhuhr-minutes arg             minutes after mid-way for calculating Dhuhr prayer time\n"
	      " ** --maghrib-minutes arg           minutes after sunset for calculating Maghrib prayer time\n"
	      " ** --isha-minutes arg              minutes after Maghrib for calculating Isha prayer time\n"
	      " ** --imsak-angle arg               angle for calculating Imsak time\n"
	      " ** --fajr-angle arg                angle for calculating Fajr prayer time\n"
	      " ** --maghrib-angle arg             angle for calculating Maghrib prayer time\n"
	      " ** --isha-angle arg                angle for calculating Isha prayer time\n"
	      "\n"
	      "  * These options are required\n"
	      " ** By providing any of these options the calculation method is set to custom\n"
	      "\n"
	      " Possible arguments for --calc-method\n"
	      "    mwl         Muslim World League\n"
	      "    isna        Islamic Society of North America\n"
	      "    jafari      Shia Ithna-Ashari, Leva Institute, Qum\n"
	      "    karachi     University of Islamic Sciences, Karachi\n"
	      "    makkah      Umm Al-Qura University, Makkah\n"
	      "    egypt       Egyptian General Authority of Survey\n"
	      "    tehran      Institute of Geophysics, University of Tehran\n"
	      "    custom      Custom Setting\n"
	      "\n"
	      " Possible arguments for --asr-juristics-method\n"
	      "    standard      Shafi`i, Maliki, Ja`fari, Hanbali\n"
	      "    hanafi        Hanafi\n"
	      "\n"
	      " Possible arguments for --high-lats-method\n"
	      "    none          No adjustment\n"
	      "    midnight      Middle of night\n"
	      "    oneseventh    1/7th of night\n"
	      "    anglebased    Angle/60th of night\n"
	      , stderr);
}              

int main(int argc, char* argv[])
{
	PrayerTimes prayer_times;
	double latitude = NAN;
	double longitude = NAN;
	double elevation = 0;
	time_t date = time(NULL);
	double timezone = NAN;

	// Parse options
	for (;;)
	{
		static option long_options[] =
		{
			{ "help",                 no_argument,       NULL, 'h' },
			{ "version",              no_argument,       NULL, 'v' },
			{ "date",                 required_argument, NULL, 'd' },
			{ "timezone",             required_argument, NULL, 'z' },
			{ "latitude",             required_argument, NULL, 'l' },
			{ "longitude",            required_argument, NULL, 'n' },
			{ "elevation",            required_argument, NULL, 'n' },
			{ "calc-method",          required_argument, NULL, 'c' },
			{ "asr-juristics-method", required_argument, NULL, 'a' },
			{ "high-lats-method",     required_argument, NULL, 'i' },
			{ "imsak-minutes",        required_argument, NULL, 0   },
			{ "dhuhr-minutes",        required_argument, NULL, 0   },
			{ "maghrib-minutes",      required_argument, NULL, 0   },
			{ "isha-minutes",         required_argument, NULL, 0   },
			{ "imsak-angle",          required_argument, NULL, 0   },
			{ "fajr-angle",           required_argument, NULL, 0   },
			{ "maghrib-angle",        required_argument, NULL, 0   },
			{ "isha-angle",           required_argument, NULL, 0   },
			{ 0, 0, 0, 0 }
		};

		enum	// long options missing a short form
		{
			IMSAK_MINUTES = 10,
			DHUHR_MINUTES,
			MAGHRIB_MINUTES,
			ISHA_MINUTES,
			IMSAK_ANGLE,
			FAJR_ANGLE,
			MAGHRIB_ANGLE,
			ISHA_ANGLE,
		};

		int option_index = 0;
		int c = getopt_long(argc, argv, "hvd:z:l:n:e:c:a:i:", long_options, &option_index);

		if (c == -1)
			break;		// Last option

		if (!optarg && c != 'h' && c != 'v')
		{
			fprintf(stderr, "Error: %s option requires an argument\n", long_options[option_index].name);
			return 2;
		}

		switch (c)
		{
			case 0:
				double arg;
				if (sscanf(optarg, "%lf", &arg) != 1)
				{
					fprintf(stderr, "Error: Invalid number '%s'\n", optarg);
					return 2;
				}
				switch (option_index)
				{
					case IMSAK_MINUTES:
						prayer_times.set_minutes(prayertimes::Imsak, arg);
						break;
					case DHUHR_MINUTES:
						prayer_times.set_minutes(prayertimes::Dhuhr, arg);
						break;
					case MAGHRIB_MINUTES:
						prayer_times.set_minutes(prayertimes::Maghrib, arg);
						break;
					case ISHA_MINUTES:
						prayer_times.set_minutes(prayertimes::Isha, arg);
						break;
					case IMSAK_ANGLE:
						prayer_times.set_angle(prayertimes::Imsak, arg);
						break;
					case FAJR_ANGLE:
						prayer_times.set_angle(prayertimes::Fajr, arg);
						break;
					case MAGHRIB_ANGLE:
						prayer_times.set_angle(prayertimes::Maghrib, arg);
						break;
					case ISHA_ANGLE:
						prayer_times.set_angle(prayertimes::Isha, arg);
						break;
					default:
						fprintf(stderr, "Error: Invalid command line option\n");
						return 2;
				}
				break;
			case 'h':		// --help
				print_help(stdout);
				return 0;
			case 'v':		// --version
				puts(PROG_NAME_FRIENDLY " " PROG_VERSION);
				return 0;
			case 'd':		// --date
			{
				tm* new_date = getdate(optarg);
				if (!new_date)
				{
					fprintf(stderr, "Error: Failed to parse '%s' as date\n", optarg);
					return 2;
				}
				date = mktime(new_date);
				break;
			}
			case 'z':		// --timezone
				if (sscanf(optarg, "%lf", &timezone) != 1)
				{
					fprintf(stderr, "Error: Invalid timezone '%s'\n", optarg);
					return 2;
				}
				break;
			case 'l':		// --latitude
				if (sscanf(optarg, "%lf", &latitude) != 1)
				{
					fprintf(stderr, "Error: Invalid latitude '%s'\n", optarg);
					return 2;
				}
				break;
			case 'n':		// --longitude
				if (sscanf(optarg, "%lf", &longitude) != 1)
				{
					fprintf(stderr, "Error: Invalid longitude '%s'\n", optarg);
					return 2;
				}
				break;
			case 'e':		// --elevation
				if (sscanf(optarg, "%lf", &elevation) != 1)
				{
					fprintf(stderr, "Error: Invalid elevation '%s'\n", optarg);
					return 2;
				}
				break;
			case 'c':		// --calc-method
				int i;
				for (i = 0; i < prayertimes::CalculationMethodsCount; ++i)
					if (strcasecmp(optarg, CalculationMethodName[i]) == 0)
					{
						prayer_times.set_calc_method(static_cast<prayertimes::CalculationMethod>(i));
						break;
					}
				if (i == prayertimes::CalculationMethodsCount)		// If none of method names has matched
				{
					fprintf(stderr, "Error: Unknown calculation method '%s'\n", optarg);
					return 2;
				}
				break;
			case 'a':		// --asr-juristics-method
				if (strcmp(optarg, "standard") == 0)
					prayer_times.settings.asr_juristics_method = prayertimes::StandardAsr;
				else if (strcmp(optarg, "hanafi") == 0)
					prayer_times.settings.asr_juristics_method = prayertimes::HanafiAsr;
				else
				{
					fprintf(stderr, "Error: Unknown Asr juristics method '%s'\n", optarg);
					return 2;
				}
				break;
			case 'i':		// --high-lats-method
				if (strcmp(optarg, "none") == 0)
					prayer_times.settings.high_latitudes_method = prayertimes::None;
				else if (strcmp(optarg, "midnight") == 0)
					prayer_times.settings.high_latitudes_method = prayertimes::NightMiddle;
				else if (strcmp(optarg, "oneseventh") == 0)
					prayer_times.settings.high_latitudes_method = prayertimes::OneSeventh;
				else if (strcmp(optarg, "anglebased") == 0)
					prayer_times.settings.high_latitudes_method = prayertimes::AngleBased;
				else
				{
					fprintf(stderr, "Error: Unknown method '%s'\n", optarg);
					return 2;
				}
				break;
			default:
				fprintf(stderr, "Error: Unknown option '%c'\n", c);
				print_help(stderr);
				return 2;
		}
	}

	if (isnan(latitude) || isnan(longitude))
	{
		fprintf(stderr, "Error: You must provide both latitude and longitude\n");
		return 2;
	}

	fputs(PROG_NAME_FRIENDLY " " PROG_VERSION "\n\n", stderr);

	if (isnan(timezone))
		timezone = PrayerTimes::get_timezone(date);

	double times[prayertimes::TimesCount];
	fprintf(stderr, "date          : %s", ctime(&date));
	fprintf(stderr, "timezone      : %.1lf\n", timezone);
	fprintf(stderr, "latitude      : %.5lf\n", latitude);
	fprintf(stderr, "longitude     : %.5lf\n", longitude);
	fprintf(stderr, "elevation     : %.5lf\n", elevation);
	fprintf(stderr, "method        : %s\n", CalculationMethodName[prayer_times.get_calc_method()]);
	putc('\n', stderr);
	prayer_times.get_prayer_times(date, latitude, longitude, elevation, timezone, times);
	for (int i = 0; i < prayertimes::TimesCount; ++i)
	{
		int seconds = times[i];
		int hours = seconds / 3600;
		seconds %= 3600;
		int minutes = seconds / 60;
		seconds %= 60;
		const char* the_next_day = "";
		if (hours >= 24)
		{
			hours -= 24;
			the_next_day = " the next day";
		}
		printf("%8s : %.2u:%.2u:%.2u%s\n", TimeName[i], hours, minutes, seconds, the_next_day);
	}
	return 0;
}
