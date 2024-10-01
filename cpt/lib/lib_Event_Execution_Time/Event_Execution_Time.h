#pragma once

#include <vector>
#include <string>
#include <chrono>


class Event
{
public:
	Event(
		std::string tag
	);

	void init();
	void start();
	void stop();
	
	//----------------------------------------------------------
	// vars
	//----------------------------------------------------------
	std::chrono::high_resolution_clock::time_point time_start, time_stop;
	std::chrono::nanoseconds duration;
	//std::chrono::microseconds duration;
	//using picoseconds = std::chrono::duration<long long, std::pico>;

	std::string tag;
	int count_occurrence;
};


class Event_Execution_Time
{
public:
	Event_Execution_Time(
		bool active
	);
	Event_Execution_Time(int num_events = 0,
		std::string default_tag = "default_tag",
		std::string output_fname = "./event_execution_time_performance",
		bool active = true
	);
	
	void start_event(int sn_event=-1, std::string tag="default_tag");
	void start_event(int sn_event = -1);
	void start_event(std::string tag);
	void stop_event(int sn_event);
	void stop_event(std::string tag);
	void write();
	int get_sn_event_from_tag(std::string tag);

	//----------------------------------------------------------
	// vars
	//----------------------------------------------------------
	std::vector<Event> event;
	std::string default_tag;
	std::string output_fname;
	bool active;
};

