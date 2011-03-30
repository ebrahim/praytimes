/*-------------------- In the name of God ----------------------*\
 
    libprayertimes 1.0
    Islamic Prayer Times Calculator library
		 
Developed by:
  Mohammad Ebrahim Mohammadi Panah <ebrahim at mohammadi dot ir>

Based on a LGPL 3.0 JavaScript code by Hamid Zarrabi-Zadeh

------------------------------------------------------------------

Copyright 2009-2011, Mohammad Ebrahim Mohammadi Panah

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License version 3 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You can get a copy of the GNU General Public License from
http://www.gnu.org/copyleft/lgpl.html

------------- Copyright block of original JS library -------------

PrayTimes.js: Prayer Times Calculator (ver 2.1)
Copyright (C) 2007-2011 PrayTimes.org

Developer: Hamid Zarrabi-Zadeh
License: GNU LGPL v3.0

TERMS OF USE:
	Permission is granted to use this code, with or 
	without modification, in any website or application 
	provided that credit is given to the original work 
	with a link back to PrayTimes.org.

This program is distributed in the hope that it will 
be useful, but WITHOUT ANY WARRANTY. 

PLEASE DO NOT REMOVE THIS COPYRIGHT BLOCK.
 
------------------------ Help and Manual -------------------------

User's Manual: 
http://praytimes.org/manual

Calculation Formulas: 
http://praytimes.org/calculation

\*--------------------------------------------------------------*/

#include <array>
#include <cstdio>
#include <cmath>
#include <string>

namespace prayertimes
{

enum
{
	VERSION_MAJOR = 1,
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
};

enum MidnightMethod
{
	StandardMidnight,   // Mid Sunset to Sunrise
	JafariMidnight,     // Mid Maghrib to Fajr
};

// Adjust methods for higher latitudes
enum HighLatMethod
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
	bool   fajr_is_minutes;
	double fajr;
	bool   dhuhr_is_minutes;
	double dhuhr;
	bool   maghrib_is_minutes;
	double maghrib;
	bool   isha_is_minutes;
	double isha;
	MidnightMethod midnight_method;
};

class PrayerTimes
{
public:
	typedef array<double, TimesCount> TimeOffsets;

	PrayerTimes(CalculationMethod calc_method = MWL,
			AsrJuristicMethod asr_juristic = StandardAsr,
			HighLatMethod adjust_high_lats = NightMiddle)
	: asr_juristic(asr_juristic)
	, adjust_high_lats(adjust_high_lats)
	{
		method_params[MWL]     = { true, 10.0, false, 18.0, true, 0, true,  0.0, false, 17.0, StandardMidnight };
		method_params[ISNA]    = { true, 10.0, false, 15.0, true, 0, true,  0.0, false, 15.0, StandardMidnight };
		method_params[Egypt]   = { true, 10.0, false, 19.5, true, 0, true,  0.0, false, 17.5, StandardMidnight };
		method_params[Makkah]  = { true, 10.0, false, 19.0, true, 0, true,  0.0, true,  90.0, StandardMidnight };
		method_params[Karachi] = { true, 10.0, false, 18.0, true, 0, true,  0.0, false, 18.0, StandardMidnight };
		method_params[Tehran]  = { true, 10.0, false, 17.7, true, 0, false, 4.5, false, 14.0, JafariMidnight   };
		method_params[Jafari]  = { true, 10.0, false, 16.0, true, 0, false, 4.0, false, 14.0, JafariMidnight   };
		method_params[Custom]  = { true, 10.0, false, 15.0, true, 0, true,  0.0, false, 15.0, StandardMidnight };

		for (int i = 0; i < TimesCount; ++i)
			times_offset[i] = 0;

		set_calc_method(calc_method);
	}

	// Get current calculation method
	CalculationMethod get_calc_method()
	{
		return calc_method;
	}

	// Set the calculation method
	void set_calc_method(CalculationMethod new_calc_method)
	{
		calc_method = new_calc_method;
		config = method_params[calc_method % CalculationMethodsCount];		// Mod is for safety in case of invalid input
	}

	// Get current time offsets
	TimeOffsets get_time_offsets()
	{
		return time_offset;
	}

	// Set time offsets
	void set_time_offsets(Times time, double new_offset)
	{
		time_offsets[time] = new_offset;
	}

	void set_time_offsets(TimeOffsets new_time_offsets)
	{
		times_offset = new_time_offsets;
	}

	// Get current setting
	MethodConfig get_config()
	{
		return config;
	}

	// return prayer times for a given date
	getTimes: function(date, coords, timezone, dst, format) {
		lat = 1* coords[0]; 
		lng = 1* coords[1];
		elv = coords[2] ? 1* coords[2] : 0;
		timeFormat = format || timeFormat;
		if (date.constructor === Date)
			date = [date.getFullYear(), date.getMonth()+ 1, date.getDate()];
		timeZone = 1* timezone+ (dst ? 1 : 0);
		jDate = this.julian(date[0], date[1], date[2])- lng/ (15* 24);
		
		return this.computeTimes();
	},


	// convert float time to the given format (see timeFormats)
	getFormattedTime: function(time, format, suffixes) {
		if (isNaN(time))
			return InvalidTime;
		if (format == 'Float') return time;
		suffixes = suffixes || timeSuffixes;

		time = DMath.fixHour(time+ 0.5/ 60);  // add 0.5 minutes to round
		var hours = Math.floor(time); 
		var minutes = Math.floor((time- hours)* 60);
		var suffix = (format == '12h') ? suffixes[hours < 12 ? 0 : 1] : '';
		hours = (format == '24h') ? this.twoDigitsFormat(hours) : ((hours+ 12 -1)% 12+ 1);
		return hours+ ':'+ this.twoDigitsFormat(minutes)+ (suffix ? ' '+ suffix : '');
	},


	//---------------------- Calculation Functions -----------------------


	// compute mid-day time
	midDay: function(time) {
		var eqt = this.sunPosition(jDate+ time).equation;
		var noon = DMath.fixHour(12- eqt);
		return noon;
	},


	// compute the time at which sun reaches a specific angle below horizon
	sunAngleTime: function(angle, time, direction) {
		var decl = this.sunPosition(jDate+ time).declination;
		var noon = this.midDay(time);
		var t = 1/15* DMath.arccos((-DMath.sin(angle)- DMath.sin(decl)* DMath.sin(lat))/ 
				(DMath.cos(decl)* DMath.cos(lat)));
		return noon+ (direction == 'ccw' ? -t : t);
	},


	// compute asr time 
	asrTime: function(factor, time) { 
		var decl = this.sunPosition(jDate+ time).declination;
		var angle = -DMath.arccot(factor+ DMath.tan(Math.abs(lat- decl)));
		return this.sunAngleTime(angle, time);
	},


	// compute declination angle of sun and equation of time
	// Ref: http://aa.usno.navy.mil/faq/docs/SunApprox.php
	sunPosition: function(jd) {
		var D = jd - 2451545.0;
		var g = DMath.fixAngle(357.529 + 0.98560028* D);
		var q = DMath.fixAngle(280.459 + 0.98564736* D);
		var L = DMath.fixAngle(q + 1.915* DMath.sin(g) + 0.020* DMath.sin(2*g));

		var R = 1.00014 - 0.01671* DMath.cos(g) - 0.00014* DMath.cos(2*g);
		var e = 23.439 - 0.00000036* D;

		var RA = DMath.arctan2(DMath.cos(e)* DMath.sin(L), DMath.cos(L))/ 15;
		var eqt = q/15 - DMath.fixHour(RA);
		var decl = DMath.arcsin(DMath.sin(e)* DMath.sin(L));

		return {declination: decl, equation: eqt};
	},


	// convert Gregorian date to Julian day
	// Ref: Astronomical Algorithms by Jean Meeus
	julian: function(year, month, day) {
		if (month <= 2) {
			year -= 1;
			month += 12;
		};
		var A = Math.floor(year/ 100);
		var B = 2- A+ Math.floor(A/ 4);

		var JD = Math.floor(365.25* (year+ 4716))+ Math.floor(30.6001* (month+ 1))+ day+ B- 1524.5;
		return JD;
	},

	
	//---------------------- Compute Prayer Times -----------------------


	// compute prayer times at given julian date
	computePrayerTimes: function(times) {
		times = this.dayPortion(times);
		var params  = setting;
		
		var imsak   = this.sunAngleTime(this.eval(params.imsak), times.imsak, 'ccw');
		var fajr    = this.sunAngleTime(this.eval(params.fajr), times.fajr, 'ccw');
		var sunrise = this.sunAngleTime(this.riseSetAngle(), times.sunrise, 'ccw');  
		var dhuhr   = this.midDay(times.dhuhr);
		var asr     = this.asrTime(this.asrFactor(params.asr), times.asr);
		var sunset  = this.sunAngleTime(this.riseSetAngle(), times.sunset);;
		var maghrib = this.sunAngleTime(this.eval(params.maghrib), times.maghrib);
		var isha    = this.sunAngleTime(this.eval(params.isha), times.isha);

		return {
			imsak: imsak, fajr: fajr, sunrise: sunrise, dhuhr: dhuhr, 
			asr: asr, sunset: sunset, maghrib: maghrib, isha: isha
		};
	},


	// compute prayer times 
	computeTimes: function() {
		// default times
		var times = { 
			imsak: 5, fajr: 5, sunrise: 6, dhuhr: 12, 
			asr: 13, sunset: 18, maghrib: 18, isha: 18
		};

		// main iterations
		for (var i=1 ; i<=numIterations ; i++) 
			times = this.computePrayerTimes(times);

		times = this.adjustTimes(times);
		
		// add midnight time
		times.midnight = (setting.midnight == 'Jafari') ? 
				times.maghrib+ this.timeDiff(times.maghrib, times.fajr)/ 2 :
				times.sunset+ this.timeDiff(times.sunset, times.sunrise)/ 2;

		times = this.tuneTimes(times);
		return this.modifyFormats(times);
	},


	// adjust times 
	adjustTimes: function(times) {
		var params = setting;
		for (var i in times)
			times[i] += timeZone- lng/ 15;
			
		if (params.highLats != 'None')
			times = this.adjustHighLats(times);
			
		if (this.isMin(params.imsak))
			times.imsak = times.fajr- this.eval(params.imsak)/ 60;
		if (this.isMin(params.maghrib))
			times.maghrib = times.sunset+ this.eval(params.maghrib)/ 60;
		if (this.isMin(params.isha))
			times.isha = times.maghrib+ this.eval(params.isha)/ 60;
		times.dhuhr += this.eval(params.dhuhr)/ 60; 

		return times;
	},


	// get asr shadow factor
	asrFactor: function(asrParam) {
		var factor = {Standard: 1, Hanafi: 2}[asrParam];
		return factor || this.eval(asrParam);
	},


	// return sun angle for sunset/sunrise
	riseSetAngle: function() {
		//var earthRad = 6371009; // in meters
		//var angle = DMath.arccos(earthRad/(earthRad+ elv));
		var angle = 0.0347* Math.sqrt(elv); // an approximation
		return 0.833+ angle;
	},


	// apply offsets to the times
	tuneTimes: function(times) {
		for (var i in times)
			times[i] += offset[i]/ 60; 
		return times;
	},


	// convert times to given time format
	modifyFormats: function(times) {
		for (var i in times)
			times[i] = this.getFormattedTime(times[i], timeFormat); 
		return times;
	},


	// adjust times for locations in higher latitudes
	adjustHighLats: function(times) {
		var params = setting;
		var nightTime = this.timeDiff(times.sunset, times.sunrise); 

		times.imsak = this.adjustHLTime(times.imsak, times.sunrise, this.eval(params.imsak), nightTime, 'ccw');
		times.fajr  = this.adjustHLTime(times.fajr, times.sunrise, this.eval(params.fajr), nightTime, 'ccw');
		times.isha  = this.adjustHLTime(times.isha, times.sunset, this.eval(params.isha), nightTime);
		times.maghrib = this.adjustHLTime(times.maghrib, times.sunset, this.eval(params.maghrib), nightTime);
		
		return times;
	},

	// adjust a time for higher latitudes
	adjustHLTime: function(time, base, angle, night, direction) {
		var portion = this.nightPortion(angle, night);
		var timeDiff = (direction == 'ccw') ? 
			this.timeDiff(time, base):
			this.timeDiff(base, time);
		if (isNaN(time) || timeDiff > portion) 
			time = base+ (direction == 'ccw' ? -portion : portion);
		return time;
	},

	// the night portion used for adjusting times in higher latitudes
	nightPortion: function(angle, night) {
		var method = setting.highLats;
		var portion = 1/2 // MidNight
		if (method == 'AngleBased')
			portion = 1/60* angle;
		if (method == 'OneSeventh')
			portion = 1/7;
		return portion* night;
	},


	// convert hours to day portions 
	dayPortion: function(times) {
		for (var i in times)
			times[i] /= 24;
		return times;
	},



	//---------------------- Misc Functions -----------------------


	// convert given string into a number
	eval: function(str) {
		return 1* (str+ '').split(/[^0-9.+-]/)[0];
	},


	// detect if input contains 'min'
	isMin: function(arg) {
		return (arg+ '').indexOf('min') != -1;
	},


	// compute the difference between two times 
	timeDiff: function(time1, time2) {
		return DMath.fixHour(time2- time1);
	},


	// add a leading 0 if necessary
	twoDigitsFormat: function(num) {
		return (num <10) ? '0'+ num : num;
	}
	
}
}



//---------------------- Degree-Based Math Class -----------------------


var DMath = {

	dtr: function(d) { return (d * Math.PI) / 180.0; },
	rtd: function(r) { return (r * 180.0) / Math.PI; },

	sin: function(d) { return Math.sin(this.dtr(d)); },
	cos: function(d) { return Math.cos(this.dtr(d)); },
	tan: function(d) { return Math.tan(this.dtr(d)); },

	arcsin: function(d) { return this.rtd(Math.asin(d)); },
	arccos: function(d) { return this.rtd(Math.acos(d)); },
	arctan: function(d) { return this.rtd(Math.atan(d)); },

	arccot: function(x) { return this.rtd(Math.atan(1/x)); },
	arctan2: function(y, x) { return this.rtd(Math.atan2(y, x)); },

	fixAngle: function(a) { return this.fix(a, 360); },
	fixHour:  function(a) { return this.fix(a, 24 ); },

	fix: function(a, b) { 
		a = a- b* (Math.floor(a/ b));
		return (a < 0) ? a+ b : a;
	}
}

/* ---------------------- Public Variables --------------------- */

	std::array<MethodConfig, CalculationMethodsCount> method_params;

private:
/* ---------------------- Private Variables -------------------- */

	MethodConfig config;
	CalculationMethod calc_method;		// caculation method
	AsrJuristicMethod asr_juristic;		// Juristic method for Asr
	HighLatMethod adjust_high_lats;		// adjusting method for higher latitudes

	TimeOffsets time_offsets;

	double latitude;
	double longitude;
	double elevation;
	double timezone;
	double julian_date;

/* --------------------- Technical Settings -------------------- */

	static const int NUM_ITERATIONS = 1;		// Number of iterations needed to compute times
};

}
