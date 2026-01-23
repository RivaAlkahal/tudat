#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

bool ensureDirectoryExists(const std::string& path)
{
    struct stat info;
    if(stat(path.c_str(), &info)!=0)
    {
	return mkdir(path.c_str(), 0755)==0;
    }
    return (info.st_mode & S_IFDIR);
}

int main ()
{
    const double initialTime =-365.0*24.0*86400.0;
    // const double finalTime =18250.0*86400.0;
	const double finalTime = -365.0*22.0*86400.0;
    // int days = 18250;
	int days = finalTime/86400.0 * -1;
    const double arcLength = 5.0;
    const double arcDuration = arcLength*86400.0;
    //const int numArcs = 3650;
	const int numArcs = days/arcLength;
	const int hoursperday = 10;

    double arcOverlap = 240.0;

    const std::string configDir = "configs_test";

    if (!ensureDirectoryExists(configDir))
    {
	std::cerr << "Could not create config directory: " << configDir << std::endl;
	return 1;
    }

    for (int i =0; i<numArcs;++i)
    {
		double arcStart = initialTime + i*arcDuration;
		double arcEnd = arcStart + arcDuration;
		std::vector<double> initialTimesObs;
		std::vector<double> finalTimesObs;
		for (int day = 0; day<arcDuration; ++day){
		    initialTimesObs.push_back(arcStart + 1800 + day *24*3600);
		    double finalTime = initialTimesObs[day]+ hoursperday*3600;
		    if (finalTime > arcEnd){
			break;
		    }
		    finalTimesObs.push_back(finalTime);
		}
		json config;
		config["arc_index"]=-1*i;
		config["arc_start"]=arcStart;
		config["arc_end"]=arcEnd;
		config["obs_start_times"] = initialTimesObs;
		config["obs_end_times"] = finalTimesObs;

		std::ofstream out(configDir + "/config_arc_test-" + std::to_string(i) + ".json");
		if (!out)
		{
		    std::cerr << "failed to write config for arc" << i << "\n";
		    continue;
		}
		out << config.dump(2);
		out.close();
    }
    std::cout<<"arc times and configs created"<<std::endl;
    return 0;
}