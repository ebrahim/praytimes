/*-------------------------- In the name of God ----------------------------*\

    libprayertimes 2.0
    Islamic prayer times calculator library
    Based on PrayTimes 2.3 JavaScript library

----------------------------- Copyright Block --------------------------------

Copyright (C) 2007-2010 PrayTimes.org

Developed By: Ebrahim Mohammadi <ebrahim at mohammadi dot ir>
Based on a JavaScript Code By: Hamid Zarrabi-Zadeh

License: GNU General Public License, ver 3

TERMS OF USE:
    Permission is granted to use this code, with or
    without modification, in any website or application
    provided that credit is given to the original work
    with a link back to PrayTimes.org.

This program is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY.

PLEASE DO NOT REMOVE THIS COPYRIGHT BLOCK.

------------------------------------------------------------------------------

User's Manual:
http://praytimes.org/manual

Calculating Formulas:
http://praytimes.org/calculation

Code Repository:
http://code.ebrahim.ir/prayertimes/

\*--------------------------------------------------------------------------*/

#include <utility>
#include <cmath>
#include <ctime>

namespace prayertimes
{

enum
{
	VERSION_MAJOR = 2,
	VERSION_MINOR = 0,
};

/* ----------------------- Constants ----------------------- */

// Times
enum Times
{
	Imsak,
	Fajr,
	Sunrise,
	Dhuhr,
	Asr,
	Sunset,
	Maghrib,
	Isha,
	Midnight,

	TimesCount
};

enum AsrJuristicsMethod
{
	StandardAsr,		// Shafi`i, Maliki, Ja`fari, Hanbali
	HanafiAsr,
	MinutesAsr,
};

enum MidnightMethod
{
	StandardMidnight,   // Mid Sunset to Sunrise
	JafariMidnight,     // Mid Sunset to Fajr
};

// Adjust methods for higher latitudes
enum HighLatitudeMethod
{
	NightMiddle,   // Middle of night
	AngleBased,    // Angle/60th of night
	OneSeventh,    // 1/7th of night
	None,          // No adjustment
};

// Calculation methods
enum CalculationMethod
{
	MWL,        // Muslim World League
	ISNA,       // Islamic Society of North America
	Egypt,      // Egyptian General Authority of Survey
	Makkah,     // Umm Al-Qura University, Makkah
	Karachi,    // University of Islamic Sciences, Karachi
	Jafari,     // Shia Ithna-Ashari, Leva Institute, Qum
	Tehran,     // Institute of Geophysics, University of Tehran
	Custom,     // Custom Setting

	CalculationMethodsCount
};

struct MethodConfig
{
	bool   imsak_is_minutes;
	double imsak;
	bool   fajr_is_minutes;		// Unused
	double fajr;
	bool   dhuhr_is_minutes;		// Unused
	double dhuhr;
	bool   maghrib_is_minutes;
	double maghrib;
	bool   isha_is_minutes;
	double isha;
	MidnightMethod midnight_method;
};

struct Settings : public MethodConfig
{
	AsrJuristicsMethod asr_juristics_method;		// Juristic method for Asr
	double asr;		// Asr minutes if method is minutes
	HighLatitudeMethod high_latitudes_method;		// adjusting method for higher latitudes

	Settings& operator=(const MethodConfig& method_config)
	{
		reinterpret_cast<MethodConfig&>(*this) = method_config;
		return *this;
	}
};

//---------------------- Degree-Based Math Class -----------------------

class DMath
{
public:
	static double dtr(double d) { return (d * M_PI) / 180.0; }
	static double rtd(double r) { return (r * 180.0) / M_PI; }

	static double sin(double d) { return ::sin(dtr(d)); }
	static double cos(double d) { return ::cos(dtr(d)); }
	static double tan(double d) { return ::tan(dtr(d)); }

	static double arcsin(double d) { return rtd(::asin(d)); }
	static double arccos(double d) { return rtd(::acos(d)); }
	static double arctan(double d) { return rtd(::atan(d)); }

	static double arccot(double x) { return rtd(::atan(1/x)); }
	static double arctan2(double y, double x) { return rtd(::atan2(y, x)); }

	static double fix_angle(double a) { return fix(a, 360.0); }
	static double fix_hour(double a) { return fix(a, 24.0); }

	static double fix(double a, double b)
	{ 
		a = a - b * (::floor(a / b));
		return (a < 0.0) ? a + b : a;
	}
};

class PrayerTimes
{
public:
	PrayerTimes(CalculationMethod calc_method = CalculationMethod::MWL,
			AsrJuristicsMethod asr_juristics_method = StandardAsr, double asr = 0.0,		// Set asr if the method is minutes
			HighLatitudeMethod high_latitudes_method = NightMiddle)
	{
		method_params[MWL]     = { true, 10.0, false, 18.0, true, 0.0, true,  0.0, false, 17.0, StandardMidnight };
		method_params[ISNA]    = { true, 10.0, false, 15.0, true, 0.0, true,  0.0, false, 15.0, StandardMidnight };
		method_params[Egypt]   = { true, 10.0, false, 19.5, true, 0.0, true,  0.0, false, 17.5, StandardMidnight };
		method_params[Makkah]  = { true, 10.0, false, 19.0, true, 0.0, true,  0.0, true,  90.0, StandardMidnight };
		method_params[Karachi] = { true, 10.0, false, 18.0, true, 0.0, true,  0.0, false, 18.0, StandardMidnight };
		method_params[Tehran]  = { true, 10.0, false, 17.7, true, 0.0, false, 4.5, false, 14.0, JafariMidnight   };
		method_params[Jafari]  = { true, 10.0, false, 16.0, true, 0.0, false, 4.0, false, 14.0, JafariMidnight   };
		method_params[Custom]  = { true, 10.0, false, 15.0, true, 0.0, true,  0.0, false, 15.0, StandardMidnight };

		for (int i = 0; i < TimesCount; ++i)
			time_offsets[i] = 0.0;

		settings.asr_juristics_method = asr_juristics_method;
		settings.asr = asr;
		settings.high_latitudes_method = high_latitudes_method;

		set_calc_method(calc_method);
	}

	// Return prayer times for a given date
	void get_prayer_times(int year, int month, int day,
			double latitude, double longitude, double elevation,
			double timezone, double times[])
	{
		this->latitude = latitude;
		this->longitude = longitude;
		this->elevation = elevation;
		this->timezone = timezone;
		julian_date = julian(year, month, day) - longitude / (double) (15 * 24);
		compute_times(times);
	}

	// Facility function to get date as a single time_t instead of separate parts
	void get_prayer_times(time_t date, double latitude, double longitude, double elevation, double timezone, double times[])
	{
		tm* t = localtime(&date);
		get_prayer_times(1900 + t->tm_year, t->tm_mon + 1, t->tm_mday, latitude, longitude, elevation, timezone, times);
	}

	//------------------ Configuration Functions -------------------

	// Get current calculation method
	CalculationMethod get_calc_method()
	{
		return calc_method;
	}

	// Set the calculation method
	void set_calc_method(CalculationMethod new_calc_method)
	{
		calc_method = new_calc_method;
		settings = method_params[calc_method];
	}

	// Get current time offsets
	void get_time_offsets(double times[])
	{
		for (int i = 0; i < TimesCount; ++i)
			times[i] = time_offsets[i];
	}

	// Set time offsets
	void set_time_offsets(double new_time_offsets[])
	{
		for (int i = 0; i < TimesCount; ++i)
			time_offsets[i] = new_time_offsets[i];
	}

	// Set a single time offset
	void set_time_offset(Times time, double new_offset)
	{
		time_offsets[time] = new_offset;
	}

	// Adjust a time using minutes
	void set_minutes(Times time, double minutes)
	{
		switch (time)
		{
			case Imsak:
				settings.imsak_is_minutes = true;
				settings.imsak = minutes;
				break;

			case Dhuhr:
				settings.dhuhr_is_minutes = true;
				settings.dhuhr = minutes;
				break;

			case Maghrib:
				settings.maghrib_is_minutes = true;
				settings.maghrib = minutes;
				break;

			case Isha:
				settings.isha_is_minutes = true;
				settings.isha = minutes;
				break;

			default:
				return;
		}
		calc_method = Custom;
	}

	// Adjust a time using angle
	void set_angle(Times time, double angle)
	{
		switch (time)
		{
			case Imsak:
				settings.imsak_is_minutes = false;
				settings.imsak = angle;
				break;

			case Fajr:
				settings.fajr_is_minutes = false;
				settings.fajr = angle;
				break;

			case Maghrib:
				settings.maghrib_is_minutes = false;
				settings.maghrib = angle;
				break;

			case Isha:
				settings.isha_is_minutes = false;
				settings.isha = angle;
				break;

			default:
				return;
		}
		calc_method = Custom;
	}

	//-------------------- Timezone Functions --------------------

	// Compute local timezone for a specific Gregorian local timestamp
	static double get_timezone(time_t local_time)
	{
		tm* tmp = localtime(&local_time);
		tmp->tm_isdst = 0;
		time_t local = mktime(tmp);
		tmp = gmtime(&local_time);
		tmp->tm_isdst = 0;
		time_t gmt = mktime(tmp);
		return (local - gmt) / 3600.0;
	}

	// Compute local timezone for a specific Gregorian date
	// Possible values for dst:
	//   <0: Automatic daylight saving
	//    0: No daylight saving
	//   >0: Yes daylight saving
	static double get_timezone(int year, int month, int day, int dst = -1)
	{
		tm date = { 0 };
		date.tm_year = year - 1900;
		date.tm_mon = month - 1;
		date.tm_mday = day;
		date.tm_isdst = dst;
		time_t local = mktime(&date);		// Seconds since midnight Jan 1, 1970
		return get_timezone(local);
	}

protected:
	//---------------------- Calculation Functions -----------------------

	// Compute mid-day time
	double mid_day(double time)
	{
		double equation = sun_position(julian_date + time).first;
		double noon = DMath::fix_hour(12.0 - equation);
		return noon;
	}

	// Compute the time at which sun reaches a specific angle below horizon
	double sun_angle_time(double angle, double time, bool direction_is_ccw = false)
	{
		double declination = sun_position(julian_date + time).second;
		double t = DMath::arccos((-DMath::sin(angle) -
					DMath::sin(declination) * DMath::sin(latitude)) /
				(DMath::cos(declination) * DMath::cos(latitude))) / 15.0;
		double noon = mid_day(time);
		return noon + (direction_is_ccw ? -t : t);
	}

	// Compute Asr time 
	double asr_time(double factor, double time)
	{ 
		double declination = sun_position(julian_date + time).second;
		double angle = -DMath::arccot(factor + DMath::tan(::fabs(latitude - declination)));
		return sun_angle_time(angle, time);
	}

	// Compute declination angle of sun and equation of time
	// Ref: http://aa.usno.navy.mil/faq/docs/SunApprox.php
	std::pair<double, double> sun_position(double jd)
	{
		double D = jd - 2451545.0;
		double g = DMath::fix_angle(357.529 + 0.98560028 * D);
		double q = DMath::fix_angle(280.459 + 0.98564736 * D);
		double L = DMath::fix_angle(q + 1.915* DMath::sin(g) + 0.020 * DMath::sin(2 * g));

		// double R = 1.00014 - 0.01671* DMath::cos(g) - 0.00014 * DMath::cos(2 * g);
		double e = 23.439 - 0.00000036 * D;

		double RA = DMath::arctan2(DMath::cos(e) * DMath::sin(L), DMath::cos(L)) / 15.0;
		double equation = q / 15.0 - DMath::fix_hour(RA);
		double declination = DMath::arcsin(DMath::sin(e) * DMath::sin(L));
		return { equation, declination };
	}

	// convert Gregorian date to Julian day
	// Ref: Astronomical Algorithms by Jean Meeus
	double julian(int year, int month, int day)
	{
		while (month <= 2)
		{
			year -= 1;
			month += 12;
		}

		double a = floor(year / 100.0);
		double b = 2 - a + floor(a / 4.0);

		return floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + b - 1524.5;
	}

	//---------------------- Compute Prayer Times -----------------------

	// Array of times must have at least TimesCount elements

	// Compute prayer times at given julian date
	void compute_prayer_times(double times[])
	{
		day_portion(times);

		times[Imsak]   = sun_angle_time(settings.imsak, times[Imsak], true);
		times[Fajr]    = sun_angle_time(settings.fajr, times[Fajr], true);
		times[Sunrise] = sun_angle_time(rise_set_angle(), times[Sunrise], true);  
		times[Dhuhr]   = mid_day(times[Dhuhr]);
		times[Asr]     = asr_time(asr_factor(settings.asr), times[Asr]);
		times[Sunset]  = sun_angle_time(rise_set_angle(), times[Sunset]);
		times[Maghrib] = sun_angle_time(settings.maghrib, times[Maghrib]);
		times[Isha]    = sun_angle_time(settings.isha, times[Isha]);
	}

	// Compute prayer times
	void compute_times(double times[])
	{
		static const double default_times[] = { 5, 5, 6, 12, 13, 18, 18, 18 };

		for (int i = 0; i < TimesCount; ++i)
			times[i] = default_times[i];

		// Main iterations
		for (int i = 1; i <= NUM_ITERATIONS; ++i) 
			compute_prayer_times(times);

		adjust_times(times);

		// Add midnight time
		if (settings.midnight_method == JafariMidnight)
			times[Midnight] = times[Sunset] + time_diff(times[Maghrib], times[Fajr]) / 2.0;
		else
			times[Midnight] = times[Sunset] + time_diff(times[Sunset], times[Sunrise]) / 2.0;

		tune_times(times);
		modify_formats(times);
	}

	void adjust_times(double times[])
	{
		for (int i = 0; i < TimesCount; ++i)
			times[i] += timezone - longitude / 15.0;

		if (settings.high_latitudes_method != None)
			adjust_high_latitudes(times);

		if (settings.imsak_is_minutes)
			times[Imsak] = times[Fajr] - settings.imsak / 60.0;
		if (settings.maghrib_is_minutes)
			times[Maghrib] = times[Sunset] + settings.maghrib / 60.0;
		if (settings.isha_is_minutes)
			times[Isha] = times[Maghrib] + settings.isha / 60.0;
		times[Dhuhr] += settings.dhuhr / 60.0; 
	}

	// Get Asr shadow factor
	double asr_factor(double asr_param)
	{
		switch (settings.asr_juristics_method)
		{
			case StandardAsr:
				return 1.0;
			case HanafiAsr:
				return 2.0;
			default:
				return asr_param;
		}
	}

	// Return sun angle for sunset/sunrise
	double rise_set_angle()
	{
		// double earth_rad = 6371009.0;		// In meters
		// double angle = DMath::arccos(earth_rad / (earth_rad + elevation));
		double angle = 0.0347 * ::sqrt(elevation);		// An approximation
		return 0.833 + angle;
	}

	// Apply offsets to the times
	void tune_times(double times[])
	{
		for (int i = 0; i < TimesCount; ++i)
			times[i] += time_offsets[i] / 60.0; 
	}

	// Convert times from hours to seconds
	void modify_formats(double times[])
	{
		for (int i = 0; i < TimesCount; ++i)
			times[i] *= 3600.0;
	}

	// adjust times for locations in higher latitudes
	void adjust_high_latitudes(double times[])
	{
		double night_time = time_diff(times[Sunset], times[Sunrise]); 

		times[Imsak]   = adjust_high_latitude_time(times[Imsak],   times[Sunrise], settings.imsak,   night_time, true);
		times[Fajr]    = adjust_high_latitude_time(times[Fajr],    times[Sunrise], settings.fajr,    night_time, true);
		times[Isha]    = adjust_high_latitude_time(times[Isha],    times[Sunset],  settings.isha,    night_time);
		times[Maghrib] = adjust_high_latitude_time(times[Maghrib], times[Sunset],  settings.maghrib, night_time);
	}

	// adjust a time for higher latitudes
	double adjust_high_latitude_time(double time, double base, double angle, double night, bool direction_is_ccw = false)
	{
		double portion = night_portion(angle, night);
		double time_diff_value = direction_is_ccw ? time_diff(time, base) : time_diff(base, time);
		if (time_diff_value > portion) 
			time = base + (direction_is_ccw ? -portion : portion);
		return time;
	}

	// the night portion used for adjusting times in higher latitudes
	double night_portion(double angle, double night)
	{
		double portion = 0.5;		// Midnight
		if (settings.high_latitudes_method == AngleBased)
			portion = angle / 60.0;
		else if (settings.high_latitudes_method == OneSeventh)
			portion = 1.0 / 7.0;
		return portion * night;
	}

	// Convert hours to day portions 
	void day_portion(double times[])
	{
		for (int i = 0; i < TimesCount; ++i)
			times[i] /= 24.0;
	}

	//---------------------- Misc Functions -----------------------

	// Compute the difference between two times 
	double time_diff(double time1, double time2)
	{
		return DMath::fix_hour(time2 - time1);
	}

public:
/* ---------------------- Public Variables --------------------- */

	MethodConfig method_params[CalculationMethodsCount];
	Settings settings;

protected:
/* ---------------------- Private Variables -------------------- */

	CalculationMethod calc_method;
	double time_offsets[TimesCount];

	// Temporary shared variables

	double latitude;
	double longitude;
	double elevation;
	double timezone;
	double julian_date;

/* --------------------- Technical Settings -------------------- */

	static const int NUM_ITERATIONS = 1;		// Number of iterations needed to compute times
};

}
