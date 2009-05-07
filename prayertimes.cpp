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

#include <ctime>
#include <cmath>
#include <cstring>
#include <unistd.h>
#include <getopt.h>

#include "prayertimes.hpp"

#define PROG_NAME "prayertimes"
#define PROG_NAME_FRIENDLY "PrayerTimes"
#define PROG_VERSION "0.2"

#if 0
// Calculation Method Names
static const char* CalculationMethodName[] =
{
	"Jafari", 	// Ithna Ashari
	"Karachi",	// University of Islamic Sciences, Karachi
	"ISNA",   	// Islamic Society of North America (ISNA)
	"MWL",    	// Muslim World League (MWL)
	"Makkah", 	// Umm al-Qura, Makkah
	"Egypt",  	// Egyptian General Authority of Survey
	"Custom", 	// Custom Setting
};

// Juristic Method Names
static const char* JuristicMethodName[] =
{
	"Shafii",    // Shafii (standard)
	"Hanafi",    // Hanafi
};

// Adjusting Method Names for Higher Latitudes
static const char* AdjustingMethodName[] =
{
	"None",      	// No adjustment
	"MidNight",  	// middle of night
	"OneSeventh",	// 1/7th of night
	"AngleBased",	// angle/60th of night
};
#endif

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

void print_help(FILE* f, option long_options[])
{
	fputs(PROG_NAME_FRIENDLY " " PROG_VERSION "\n\n", stderr);
	fputs("Usage: " PROG_NAME " [optionals...] requireds...\n\n", stderr);
	for (option* opt = long_options; opt->name; ++opt)
	{
		fprintf(f, "    -%c --%s %s\n", opt->val ? opt->val : '-', opt->name, opt->has_arg ? "arg" : "");
		// TODO: Add description of commands and their possible arguments
	}
}

int main(int argc, char* argv[])
{
	PrayerTimes prayer_times;
	double latitude = NAN;		// 35.7061
	double longitude = NAN;		// 51.4358
	time_t date = time(NULL);
	double timezone = NAN;

	// Parse options
	for (;;)
	{
		static option long_options[] =
		{
			{ "help",                no_argument,       NULL, 'h' },
			{ "version",             no_argument,       NULL, 'v' },
			{ "date",                required_argument, NULL, 'd' },
			{ "timezone",            required_argument, NULL, 'z' },
			{ "latitude",            required_argument, NULL, 'l' },
			{ "longitude",           required_argument, NULL, 'n' },
			{ "calc-method",         required_argument, NULL, 'c' },
			{ "asr-juristic-method", required_argument, NULL, 'a' },
			{ "high-lats-method",    required_argument, NULL, 'i' },
			{ "dhuhr-minutes",       required_argument, NULL, 0   },
			{ "maghrib-minutes",     required_argument, NULL, 0   },
			{ "isha-minutes",        required_argument, NULL, 0   },
			{ "fajr-angle",          required_argument, NULL, 0   },
			{ "maghrib-angle",       required_argument, NULL, 0   },
			{ "isha-angle",          required_argument, NULL, 0   },
			{ 0, 0, 0, 0 }
		};

		enum	// long options missing a short form
		{
			DHUHR_MINUTES = 9,
			MAGHRIB_MINUTES,
			ISHA_MINUTES,
			FAJR_ANGLE,
			MAGHRIB_ANGLE,
			ISHA_ANGLE,
		};

		int option_index = 0;
		int c = getopt_long(argc, argv, "hvd:z:l:n:c:a:i:", long_options, &option_index);

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
					case DHUHR_MINUTES:
						prayer_times.set_dhuhr_minutes(arg);
						break;
					case MAGHRIB_MINUTES:
						prayer_times.set_maghrib_minutes(arg);
						break;
					case ISHA_MINUTES:
						prayer_times.set_isha_minutes(arg);
						break;
					case FAJR_ANGLE:
						prayer_times.set_fajr_angle(arg);
						break;
					case MAGHRIB_ANGLE:
						prayer_times.set_maghrib_angle(arg);
						break;
					case ISHA_ANGLE:
						prayer_times.set_isha_angle(arg);
						break;
					default:
						fprintf(stderr, "Error: Invalid command line option\n");
						return 2;
				}
				break;
			case 'h':		// --help
				print_help(stdout, long_options);
				return 0;
			case 'v':		// --version
				puts(PROG_NAME_FRIENDLY " " PROG_VERSION);
				return 0;
			case 'd':		// --date
			{
				tm* new_date = getdate(optarg);
				if (!new_date)
				{
					fprintf(stderr, "Error: Failed to parse '%s' as date (%m)\n", optarg);
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
			case 'c':		// --calc-method
				if (strcmp(optarg, "jafari") == 0)
					prayer_times.set_calc_method(PrayerTimes::Jafari);
				else if (strcmp(optarg, "karachi") == 0)
					prayer_times.set_calc_method(PrayerTimes::Karachi);
				else if (strcmp(optarg, "isna") == 0)
					prayer_times.set_calc_method(PrayerTimes::ISNA);
				else if (strcmp(optarg, "mwl") == 0)
					prayer_times.set_calc_method(PrayerTimes::MWL);
				else if (strcmp(optarg, "makkah") == 0)
					prayer_times.set_calc_method(PrayerTimes::Makkah);
				else if (strcmp(optarg, "egypt") == 0)
					prayer_times.set_calc_method(PrayerTimes::Egypt);
				else if (strcmp(optarg, "custom") == 0)
					prayer_times.set_calc_method(PrayerTimes::Custom);
				else
				{
					fprintf(stderr, "Error: Unknown method '%s'\n", optarg);
					return 2;
				}
				break;
			case 'a':		// --asr-juristic-method
				if (strcmp(optarg, "shafii") == 0)
					prayer_times.set_asr_method(PrayerTimes::Shafii);
				else if (strcmp(optarg, "hanafi") == 0)
					prayer_times.set_asr_method(PrayerTimes::Hanafi);
				else
				{
					fprintf(stderr, "Error: Unknown method '%s'\n", optarg);
					return 2;
				}
				break;
			case 'i':		// --high-lats-method
				if (strcmp(optarg, "none") == 0)
					prayer_times.set_high_lats_adjust_method(PrayerTimes::None);
				else if (strcmp(optarg, "midnight") == 0)
					prayer_times.set_high_lats_adjust_method(PrayerTimes::MidNight);
				else if (strcmp(optarg, "oneseventh") == 0)
					prayer_times.set_high_lats_adjust_method(PrayerTimes::OneSeventh);
				else if (strcmp(optarg, "anglebased") == 0)
					prayer_times.set_high_lats_adjust_method(PrayerTimes::AngleBased);
				else
				{
					fprintf(stderr, "Error: Unknown method '%s'\n", optarg);
					return 2;
				}
				break;
			default:
				fprintf(stderr, "Error: Unknown option '%c'\n", c);
				print_help(stderr, long_options);
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
		timezone = PrayerTimes::get_effective_timezone(date);

	double times[PrayerTimes::TimesCount];
	fprintf(stderr, "date          : %s", ctime(&date));
	fprintf(stderr, "timezone      : %.1lf\n", timezone);
	fprintf(stderr, "latitude      : %.5lf\n", latitude);
	fprintf(stderr, "longitude     : %.5lf\n", longitude);
	puts("");
	prayer_times.get_prayer_times(date, latitude, longitude, timezone, times);
	for (int i = 0; i < PrayerTimes::TimesCount; ++i)
		printf("%8s : %s\n", TimeName[i], PrayerTimes::float_time_to_time24(times[i]).c_str());
	return 0;
}
