
#include "Event_Execution_Time.h"
#include <iostream>
#include <fstream>
#include <string>

Event::Event(
	std::string tag
) : tag(tag)
{
	init();
}

void Event::init()
{
	auto this_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(this_time - this_time);
	count_occurrence = 0;
}

void Event::start()
{
	time_start = std::chrono::high_resolution_clock::now();
}

void Event::stop()
{
	time_stop = std::chrono::high_resolution_clock::now();
	duration += std::chrono::duration_cast<std::chrono::nanoseconds>(time_stop - time_start);
	count_occurrence++;
}




Event_Execution_Time::Event_Execution_Time(
	bool active
) :
	active(active)
{
}

Event_Execution_Time::Event_Execution_Time(
	int num_events,
	std::string default_tag,
	std::string output_fname,
	bool active
	) : 
	default_tag(default_tag), 
	output_fname(output_fname),
	active(active)
{
	if (num_events > 0)
	{
		for (int sn_event =0; sn_event <num_events; sn_event++)
			event.push_back(Event(default_tag));
	}
}

void Event_Execution_Time::start_event(int sn_event, std::string tag)
{
	if (!active)
		return;

	if (sn_event < 0)
	{
		event.push_back(Event(tag));
		event.back().start();
	}
	else
		event[sn_event].start();
}

void Event_Execution_Time::start_event(int sn_event)
{
	if (!active)
		return;

	if (sn_event < 0)
	{
		event.push_back(Event(default_tag));
		event.back().start();
	}
	else
		event[sn_event].start();
}


void Event_Execution_Time::start_event(std::string tag)
{
	if (!active)
		return;

	int sn_event = get_sn_event_from_tag(tag);
	if (sn_event < 0)
	{
		event.push_back(Event(tag));
		event.back().start();
	}
	else
		event[sn_event].start();
}



void Event_Execution_Time::stop_event(int sn_event)
{
	if (!active)
		return;

	event[sn_event].stop();
}

void Event_Execution_Time::stop_event(std::string tag)
{
	if (!active)
		return;

	int sn_event = get_sn_event_from_tag(tag);
	if (sn_event >= 0)
		event[sn_event].stop();
	else
		std::cout << "The following provided event tag is not valid: " << tag << std::endl;
}

int Event_Execution_Time::get_sn_event_from_tag(std::string tag)
{
	for (int sn_event = 0; sn_event < event.size(); sn_event++)
		if (event[sn_event].tag == tag)
			return sn_event;
	return -1;
}

void Event_Execution_Time::write()
{
	std::ofstream fout(output_fname);
	fout << "sn\ttag\tcount of occurrence\tduration\tduration_per_occurrence\n";
	for (int sn_event = 0; sn_event < event.size(); sn_event++)
	{
		fout << sn_event << "\t"
			<< event[sn_event].tag << "\t"
			<< event[sn_event].count_occurrence << "\t"
			<< event[sn_event].duration.count() << "\t";

		if (event[sn_event].count_occurrence > 0)
			fout << event[sn_event].duration.count() / event[sn_event].count_occurrence;
		else
			fout << "NA";

		fout << "\n";
	}
	fout.close();
}






