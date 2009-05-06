/*

---------- Copyright block of original JS library ------------

PrayTime: The Prayer Times Calculator (ver 1.1)
Copyright (C) 2007, Hamid Zarrabi-Zadeh

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

---------------------- Help and Manual -----------------------

PrayTime User's Manual:
http://www.cs.uwaterloo.ca/~hzarrabi/praytime/doc/manual

Calculating Formulas:
http://www.cs.uwaterloo.ca/~hzarrabi/praytime/doc/calculation

*/

#include <cmath>
#include <string>

/* --------------------- PrayTime Class ----------------------- */

class PrayerTime
{
public:
/* --------------------- User Interface ----------------------- */
/*
   get_prayer_times(date, latitude, longitude, timezone)
   get_date_prayer_times(year, month, day, latitude, longitude, timezone)

   set_calc_method(method_id)
   set_asr_method(method_id)

   set_fajr_angle(angle)
   set_maghrib_angle(angle)
   set_isha_angle(angle)
   set_dhuhr_minutes(minutes)		// minutes after mid-day
   set_maghrib_minutes(minutes)		// minutes after sunset
   set_isha_minutes(minutes)		// minutes after maghrib

   set_hl_adjust_method(method_id)		// adjust method for higher latitudes

   set_time_format(time_format)
   float_to_time_24(time)
   float_to_time_12(time)
   float_to_time_12_ns(time)
*/

	// Calculation Methods
	enum CalculationMethod
	{
		Jafari, 	// Ithna Ashari
		Karachi,	// University of Islamic Sciences, Karachi
		ISNA,   	// Islamic Society of North America (ISNA)
		MWL,    	// Muslim World League (MWL)
		Makkah, 	// Umm al-Qura, Makkah
		Egypt,  	// Egyptian General Authority of Survey
		Custom, 	// Custom Setting

		CalculationMethodsCount
	};

	// Juristic Methods
	enum JuristicMethod
	{
		Shafii,    // Shafii (standard)
		Hanafi,    // Hanafi
	};

	// Adjusting Methods for Higher Latitudes
	enum AdjustingMethod
	{
		None,      	// No adjustment
		MidNight,  	// middle of night
		OneSeventh,	// 1/7th of night
		AngleBased,	// angle/60th of night
	};

	// Time Formats
	enum TimeFormat
	{
		Time24     = 0,    // 24-hour format
		Time12     = 1,    // 12-hour format
		Time12NS   = 2,    // 12-hour format with no suffix
		Float      = 3,    // floating point number
	};

	// Time Names
	enum TimeID
	{
		Fajr,
		Sunrise,
		Dhuhr,
		Asr,
		Sunset,
		Maghrib,
		Isha,

		TimeNameSize
	};

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

	static const std::string InvalidTime = "-----";	 // The string used for invalid times

/* ------------------- Calc Method Parameters -------------------- */

	struct MethodConfig
	{
		double fajr_angle;
		bool maghrib_is_minutes;
		double maghrib_value;		// angle or minutes
		bool isha_is_minutes;
		double isha_value;		// angle or minutes
	};

	MethodConfig method_params[] =
	{
		{ 16,   false, 4, false, 14   },	// Jafari
		{ 18,   true,  0, false, 18   },	// Karachi
		{ 15,   true,  0, false, 15   },	// ISNA
		{ 18,   true,  0, false, 17   },	// MWL
		{ 19,   true,  0, true,  90   },	// Makkah
		{ 19.5, true,  0, false, 17.5 },	// Egypt
		{ 18,   true,  0, false, 17   },	// Custom
	};

/* -------------------- Interface Functions -------------------- */

	/* return prayer times for a given date */
	/*TODO*/ get_date_prayer_times(int year, int month, int day, double _latitude, double _longitude, double _timezone, double times[])
	{
		latitude = _latitude;
		longitude = _longitude;
		julian_date = get_julian_date(year, month, day) - longitude / (15 * 24);
		compute_day_times(times);
	}

#if 0
	/* return prayer times for a given date */
	get_prayer_times(date, latitude, longitude, timezone)
	{
		return get_date_prayer_times(date.get_full_year(), date.get_month() + 1, date.get_date(),
					latitude, longitude, timezone);
	}
#endif

	/* set the calculation method  */
	void set_calc_method(CalculationMethod method_id)
	{
		calc_method = method_id;
	}

	/* set the juristic method for Asr */
	void set_asr_method(JuristicMethod method_id)
	{
		asr_juristic = method_id;
	}

	/* set the angle for calculating Fajr */
	void set_fajr_angle(double angle)
	{
		method_params[Custom].fajr_angle = angle;
		calc_method = Custom;
	}

	/* set the angle for calculating Maghrib */
	void set_maghrib_angle(angle)
	{
		method_params[Custom].maghrib_is_minutes = false;
		method_params[Custom].maghrib_value = angle;
		calc_method = Custom;
	}

	/* set the angle for calculating Isha */
	void set_isha_angle(angle)
	{
		method_params[Custom].isha_is_minutes = false;
		method_params[Custom].isha_value = angle;
		calc_method = Custom;
	}

	/* set the minutes after mid-day for calculating Dhuhr */
	void set_dhuhr_minutes(minutes)
	{
		dhuhr_minutes = minutes;
	}

	/* set the minutes after Sunset for calculating Maghrib */
	void set_maghrib_minutes(minutes)
	{
		method_params[Custom].maghrib_is_minutes = true;
		method_params[Custom].maghrib_value = minutes;
		calc_method = Custom;
	}

	/* set the minutes after Maghrib for calculating Isha */
	void set_isha_minutes(minutes)
	{
		method_params[Custom].isha_is_minutes = true;
		method_params[Custom].isha_value = minutes;
		calc_method = Custom;
	}

	/* set adjusting method for higher latitudes */
	void set_hl_adjust_method(AdjustingMethod method_id)
	{
		adjust_high_lats = method_id;
	}

	/* set the time format */
	void set_time_format(TimeFormat _time_format)
	{
		time_format = _time_format;
	}
#if 0
	/* convert float hours to 24h format */
	std::string float_to_time24(double time)
	{
		if (is_nan(time))
			return InvalidTime;
		time = fix_hour(time + 0.5 / 60);  // add 0.5 minutes to round
		int hours = Math.floor(time);
		int minutes = Math.floor((time- hours)* 60);
		return this.two_digits_format(hours)+':'+ this.two_digits_format(minutes);
	}

	/* convert float hours to 12h format */
	std::string float_to_time12(time, no_suffix)
	{
		if (is_nan(time))
			return this.InvalidTime;
		time = this.fix_hour(time+ 0.5/ 60);  // add 0.5 minutes to round
		var hours = Math.floor(time);
		var minutes = Math.floor((time- hours)* 60);
		var suffix = hours >= 12 ? ' pm' : ' am';
		hours = (hours+ 12 -1)% 12+ 1;
		return hours+':'+ this.two_digits_format(minutes)+ (no_suffix ? '' : suffix);
	}

	/* convert float hours to 12h format with no suffix */
	std::string float_to_time12ns(time)
	{
		return this.float_to_time12(time, true);
	}
#endif

/* ---------------------- Time-Zone Functions ----------------------- */

	/* compute local time-zone for a specific date */
	double get_timezone(time_t local_time)
	{
		time_t gmt_time = mktime(gmtime(&local_time));
		return (local_time - gmt_time) / (double) (60 * 60);
	}

	/* return effective timezone for a given date */
	double effective_timezone(int year, int month, int day)
	{
		tm t = { 0 };
		t.tm_mday = day;
		t.tm_mon = month - 1;
		t.tm_year = year;
		return get_timezone(mktime(&t));
	}

#if 0
	/* compute base time-zone of the system */
	get_base_timezone()
	{
		return this.get_timezone(new Date(2000, 0, 15))
	}

	/* detect daylight saving in a given date */
	detect_daylight_saving(date)
	{
		return this.get_timezone(date) != this.get_base_timezone();
	}
#endif

private:
/* ---------------------- Calculation Functions ----------------------- */

	/* References: */
	/* http://www.ummah.net/astronomy/saltime   */
	/* http://aa.usno.navy.mil/faq/docs/SunApprox.html */

	typedef std::pair<double, double> DoublePair;

	/* compute declination angle of sun and equation of time */
	DoublePair sun_position(double jd)
	{
		double d = jd - 2451545.0;
		double g = fix_angle(357.529 + 0.98560028 * d);
		double q = fix_angle(280.459 + 0.98564736 * d);
		double l = fix_angle(q + 1.915 * dsin(g) + 0.020 * dsin(2 * g));

		double r = 1.00014 - 0.01671 * dcos(g) - 0.00014 * dcos(2 * g);
		double e = 23.439 - 0.00000036 * d;

		double d = darcsin(dsin(e) * dsin(_l));
		double ra = darctan2(dcos(e) * dsin(_l), dcos(l)) / 15.0;
		ra = fix_hour(ra);
		double eq_t = q / 15.0 - ra;

		return DoublePair(d, eq_t);
	}

	/* compute equation of time */
	double equation_of_time(double jd)
	{
		return sun_position(jd).second;
	}

	/* compute declination angle of sun */
	double sun_declination(double jd)
	{
		return sun_position(jd).first;
	}

	/* compute mid-day (Dhuhr, Zawal) time */
	double compute_mid_day(double t)
	{
		double t = equation_of_time(julian_date + t);
		double z = fix_hour(12 - t);
		return z;
	}

	/* compute time for a given angle G */
	double compute_time(double g, double t)
	{
		var d = sun_declination(julian_date + t);
		var z = compute_mid_day(t);
		var v = 1.0 / 15.0 * darccos((-dsin(g) - dsin(d) * dsin(latitude)) / (dcos(d) * dcos(latitude)));
		return z + (g > 90.0 ? - v :  v);
	}

	/* compute the time of Asr */
	double compute_asr(int step, double t)  // Shafii: step=1, Hanafi: step=2
	{
		var d = sun_declination(julian_date + t);
		var g = -darccot(step + dtan(fabs(latitude - d)));
		return compute_time(g, t);
	}

/* ---------------------- Compute Prayer Times ----------------------- */

	// array parameters must be at least of size TimeNameSize

	/* compute prayer times at given julian date */
	void compute_times(double times[])
	{
		day_portion(times);

		times[Fajr]    = compute_time(180.0 - method_params[calc_method].fajr_angle, times[Fajr]);
		times[Sunrise] = compute_time(180.0 - 0.833, times[Sunrise]);
		times[Dhuhr]   = compute_mid_day(times[Dhuhr]);
		times[Asr]     = compute_asr(1 + asr_juristic, times[Asr]);
		times[Sunset]  = compute_time(0.833, times[Sunset]);;
		times[Maghrib] = compute_time(method_params[calc_method].maghrib_value, times[Maghrib]);
		times[Isha]    = compute_time(method_params[calc_method].isha_value, times[Isha]);
	}


	/* compute prayer times at given julian date */
	void compute_day_times(double times[])
	{
		times = { 5, 6, 12, 13, 18, 18, 18 };		// default times

		for (int i = 0; i < num_iterations; ++i)
			compute_times(times);

		adjust_times(times);
		adjust_times_format(times);
	}


	/* adjust times in a prayer time array */
	void adjust_times(double times[])
	{
		for (int i = 0 ; i < TimeNameSize; ++i)
			times[i] += timezone - longitude / 15.0;
		times[Dhuhr] += dhuhr_minutes / 60.0;		// Dhuhr
		if (method_params[calc_method].maghrib_is_minutes)		// Maghrib
			times[Maghrib] = times[Sunset] + method_params[calc_method].maghrib_value / 60.0;
		if (method_params[calc_method].isha_is_minutes)		// Isha
			times[Isha] = times[Maghrib] + method_params[calc_method].isha_value / 60.0;

		if (adjust_high_lats != None)
			adjust_high_lat_times(times);
	}

#if 0
	/* convert times array to given time format */
	adjust_times_format(double times[])
	{
		if (time_format == Float)
			return times;
		for (var i=0; i<TimeNameSize; i++)
			if (time_format == Time12)
				times[i] = float_to_time12(times[i]);
			else if (time_format == Time12NS)
				times[i] = float_to_time12(times[i], true);
			else
				times[i] = float_to_time24(times[i]);
		return times;
	}
#endif

	/* adjust Fajr, Isha and Maghrib for locations in higher latitudes */
	void adjust_high_lat_times(double times[])
	{
		double night_time = time_diff(times[Sunset], times[Sunrise]);		// sunset to sunrise

		// Adjust Fajr
		double fajr_diff = night_portion(method_params[calc_method].fajr_angle) * night_time;
		if (is_nan(times[Fajr]) || time_diff(times[Fahr], times[Sunrise]) > fajr_diff)
			times[Fajr] = times[Sunrise] - fajr_diff;

		// Adjust Isha
		double isha_angle = method_params[calc_method].isha_is_minutes ? 18.0 : method_params[calc_method].isha_value;
		double isha_diff = night_portion(isha_angle) * night_time;
		if (is_nan(times[Isha]) || time_diff(times[Sunset], times[Isha]) > isha_diff)
			times[Isha] = times[Sunset] + isha_diff;

		// Adjust Maghrib
		double maghrib_angle = method_params[calc_method].maghrib_is_minutes ? 4.0 : method_params[calc_method].maghrib_value;
		double maghrib_diff = night_portion(maghrib_angle) * night_time;
		if (is_nan(times[Maghrib]) || time_diff(times[Sunset], times[Maghrib]) > maghrib_diff)
			times[Maghrib] = times[Sunset] + maghrib_diff;

		return times;
	}


	/* the night portion used for adjusting times in higher latitudes */
	double night_portion(double angle)
	{
		switch (adjust_high_lats)
		{
			case AngleBased:
				return 1.0 / 60.0 * angle;
			case MidNight:
				return 1.0 / 2.0;
			case OneSeventh:
				return 1.0 / 7.0;
			default:
				// Just to return something!
				// In original library nothing was returned
				// Maybe I should throw an exception
				// It must be impossible to reach here
				return 0;
		}
	}

	/* convert hours to day portions  */
	void day_portion(double times[])
	{
		for (int i = 0; i < TimeNameSize; ++i)
			times[i] /= 24.0;
	}

/* ---------------------- Misc Functions ----------------------- */

	/* compute the difference between two times  */
	double time_diff(double time1, double time2)
	{
		return fix_hour(time2 - time1);
	}

#if 0
	/* add a leading 0 if necessary */
	two_digits_format(num)
	{
		return (num <10) ? '0'+ num : num;
	}
#endif

/* ---------------------- Julian Date Functions ----------------------- */

	/* calculate julian date from a calendar date */
	double julian_date(int year, int month, int day)
	{
		if (month <= 2)
		{
			year -= 1;
			month += 12;
		}

		double a = floor(year / 100.0);
		double b = 2 - a + floor(a / 4.0);

		return floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + b - 1524.5;
	}

	/* convert a calendar date to julian date (second method) */
	double calc_julian_date(int year, int month, int day)
	{
		double j1970 = 2440588.0;
		tm date = { 0 };
		date.tm_year = year;
		date.tm_mon = month - 1;
		date.tm_dday = day;
		time_t ms = mktime(&date);		// seconds since midnight Jan 1, 1970
		double days = floor(ms / (double) (60 * 60 * 24));
		return j1970 + days - 0.5;
	}

/* ---------------------- Trigonometric Functions ----------------------- */

	/* degree sin */
	double dsin(double d)
	{
		return sin(deg2rad(d));
	}

	/* degree cos */
	double dcos(double d)
	{
		return cos(deg2rad(d));
	}

	/* degree tan */
	double dtan(double d)
	{
		return tan(deg2rad(d));
	}

	/* degree arcsin */
	double darcsin(double x)
	{
		return rad2deg(asin(x));
	}

	/* degree arccos */
	double darccos(double x)
	{
		return rad2deg(acos(x));
	}

	/* degree arctan */
	double darctan(double x)
	{
		return rad2deg(atan(x));
	}

	/* degree arctan2 */
	double darctan2(double y, double x)
	{
		return rad2deg(atan2(y, x));
	}

	/* degree arccot */
	double darccot(double x)
	{
		return rad2deg(atan(1.0 / x));
	}

	/* degree to radian */
	double deg2rad(double d)
	{
		return d * M_PI / 180.0;
	}

	/* radian to degree */
	double rad2deg(double r)
	{
		return r * 180.0 / M_PI;
	}

	/* range reduce angle in degrees. */
	double fix_angle(double a)
	{
		a = a - 360.0 * floor(a / 360.0);
		a = a < 0.0 ? a + 360.0 : a;
		return a;
	}

	/* range reduce hours to 0..23 */
	double fix_hour(double a)
	{
		a = a - 24.0 * floor(a / 24.0);
		a = a < 0.0 ? a + 24.0 : a;
		return a;
	}

private:
/* ---------------------- Private Variables -------------------- */

	CalculationMethod calc_method = Jafari;		// caculation method
	JuristicMethod asr_juristic = Shafii;		// Juristic method for Asr
	AdjustingMethod adjust_high_lats = MidNight;	// adjusting method for higher latitudes
	double dhuhr_minutes = 0;		// minutes after mid-day for Dhuhr
	TimeFormat time_format = Time24;
	double latitude;
	double longitude;
	double timezone;
	double julian_date;

/* --------------------- Technical Settings -------------------- */

	static const num_iterations = 1;		// number of iterations needed to compute times
};
