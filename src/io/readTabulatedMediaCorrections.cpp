/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readTabulatedMediaCorrections.h"

#include <map>
#include <fstream>
#include <string>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace input_output
{

AtmosphericCorrectionCspCommand::AtmosphericCorrectionCspCommand( std::vector< std::string >& cspCommand ):
    sourceSpecifier_( "" ),
    sourceId_( 0 )
{
    for ( unsigned int i = 0; i < cspCommand.size(); ++i )
    {
        if ( cspCommand.at( i ) == "MODEL" )
        {
            modelIdentifier_ = cspCommand.at( i + 1 );
            ++i;
        }
        else if ( cspCommand.at( i ) == "FROM" )
        {
            startTime_ = convertTime( cspCommand.at( i + 1 ), cspCommand.at( i + 2 ) );
            i = i + 2;
        }
        else if ( cspCommand.at( i ) == "TO" )
        {
            endTime_ = convertTime( cspCommand.at( i + 1 ), cspCommand.at( i + 2 ) );
            i = i + 2;
        }
        else if ( !boost::algorithm::find_first( cspCommand.at( i ), "BY " ).empty( )  )
        {
            computationSpecifier_ = cspCommand.at( i );
            for ( ++i; i < cspCommand.size(); ++i )
            {
                try
                {
                    computationCoefficients_.push_back( std::stod( cspCommand.at( i ) ) );
                }
                catch( ... )
                {
                    --i;
                    break;
                }
            }
        }
        else if ( cspCommand.at( i ) == "DSN" )
        {
            groundStationsId_ = cspCommand.at( i + 1 );
            ++i;
        }
        else if ( cspCommand.at( i ) == "SCID" || cspCommand.at( i ) == "QUASAR" )
        {
            sourceSpecifier_ = cspCommand.at( i );
            sourceId_ = std::stoi( cspCommand.at( i + 1 ) );
            i = i + 2;
        }
        else if ( cspCommand.at( i ) == "ADJUST" )
        {
            dataTypesIdentifier_ = cspCommand.at( i + 1 );
            ++i;
        }
        else
        {
            throw std::runtime_error( "Invalid CSP command for atmospheric correction file." );
        }

    }

}

double AtmosphericCorrectionCspCommand::convertTime( std::string yearMonthDay, std::string hoursMinutesSeconds )
{
    std::vector< std::string > vectorOfIndividualStrings;

    // Extract years, months, days
    boost::algorithm::split( vectorOfIndividualStrings, yearMonthDay, boost::algorithm::is_any_of( "/" ),
                             boost::algorithm::token_compress_on );
    if ( vectorOfIndividualStrings.size( ) != 3 )
    {
        throw std::runtime_error("Error when extracting year/month/day from CSP command: number of extracted elements (" +
            std::to_string( vectorOfIndividualStrings.size( ) ) + ") is invalid.");
    }
    int year = std::stoi( vectorOfIndividualStrings.at( 0 ) );
    if ( year >= 00 && year <= 68 )
    {
        year += 2000;
    }
    else if ( year >= 69 && year <= 99 )
    {
        year += 1900;
    }
    else
    {
        throw std::runtime_error("Error when extracting year from CSP command: " + std::to_string( year ) +
            " is not a valid year (valid years are 00 through 99).");
    }
    int month = std::stoi( vectorOfIndividualStrings.at( 1 ) );
    int day = std::stoi( vectorOfIndividualStrings.at( 2 ) );

    // Extract hours, minutes and seconds
    boost::algorithm::split( vectorOfIndividualStrings, hoursMinutesSeconds, boost::algorithm::is_any_of( ":" ),
                             boost::algorithm::token_compress_on );
    if ( vectorOfIndividualStrings.size( ) != 2 && vectorOfIndividualStrings.size( ) != 3 )
    {
        throw std::runtime_error("Error when extracting hours:minutes:seconds from CSP command: number of extracted elements (" +
            std::to_string( vectorOfIndividualStrings.size( ) ) + ") is invalid.");
    }
    int hours = std::stoi( vectorOfIndividualStrings.at( 0 ) );
    int minutes = std::stoi( vectorOfIndividualStrings.at( 1 ) );
    double seconds = 0.0;
    if ( vectorOfIndividualStrings.size( ) == 3 )
    {
        seconds = std::stod( vectorOfIndividualStrings.at( 2 ) );
    }

    return basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            year, month, day, hours, minutes, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;
}

void CspRawFile::parseCspCommands( const std::vector< std::string >& cspCommandsVector )
{

    for ( unsigned int i = 0; i < cspCommandsVector.size( ); ++i )
    {
        std::vector< std::string > cspVectorOfIndividualStrings;
        boost::algorithm::split( cspVectorOfIndividualStrings, cspCommandsVector.at( i ),
                                 boost::algorithm::is_any_of( ",()" ), boost::algorithm::token_compress_on );
        // Check if last string is empty and remove it if so
        if ( cspVectorOfIndividualStrings.back( ).empty( ) )
        {
            cspVectorOfIndividualStrings.pop_back( );
        }

        std::shared_ptr< CspCommand > cspCommand;

        // Loop over strings and find type of model
        for ( unsigned int j = 0; j < cspVectorOfIndividualStrings.size( ); ++ j )
        {
            boost::algorithm::trim( cspVectorOfIndividualStrings.at( j ) );
            if ( cspVectorOfIndividualStrings.at( j ) == "MODEL" )
            {
                if ( cspVectorOfIndividualStrings.at( j + 1 ) == "DRY NUPART" ||
                    cspVectorOfIndividualStrings.at( j + 1 ) == "WET NUPART" ||
                    cspVectorOfIndividualStrings.at( j + 1 ) == "CHPART" )
                {
                    cspCommand = std::make_shared< AtmosphericCorrectionCspCommand >( cspVectorOfIndividualStrings );
                }
                else
                {
                    throw std::runtime_error( "Error when reading CSP file: invalid model type." );
                }

                cspCommands_.push_back( cspCommand );
                break;
            }
        }
    }
}

std::vector< std::string > CspRawFile::readCspCommandsFile( const std::string& file )
{

    // Open file
    std::ifstream stream( file, std::ios_base::binary );
    if ( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening file: " + file );
    }

    // Line based parsing
    std::string line;
    std::vector< std::string > commandVector;
    bool newCommand = true;
    while ( stream.good( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Skip empty and comment lines
        if ( line.size( ) <= 0 || ( line.at( 0 ) == '#' ) )
        {
            continue;
        }

        // Search and remove in-line comments
        auto iteratorRange = boost::algorithm::find_first( line, "#" );
        if( !iteratorRange.empty( ) )
        {
            line.erase( iteratorRange.begin( ), line.end( ) );
            // Re-trim string
            boost::algorithm::trim( line );
        }

        // If continuation of previous command, merge it with previous portions
        if ( !newCommand )
        {
            // Check if it's the last part of the command
            auto iteratorRange = boost::algorithm::find_first( line, ")." );
            if ( !iteratorRange.empty( ) )
            {
                // Remove command terminator (i.e. the period)
                line.erase( ++iteratorRange.begin( ), line.end( ) );
                newCommand = true;
            }

            // Merge with previous part of command
            commandVector.back( ) += line;
        }
        // If new command, push it to vector
        else
        {
            commandVector.push_back( line );
            newCommand = false;
        }
    }

    // Close file
    stream.close();

    return commandVector;
}

std::shared_ptr< observation_models::TabulatedMediaReferenceCorrection > createReferenceCorrection(
        const double startTime,
        const double endTime,
        const std::vector< double >& coefficients,
        const std::string& computationSpecifier )
{

    std::shared_ptr< observation_models::TabulatedMediaReferenceCorrection > mediaCorrection;

    if ( computationSpecifier == "BY CONST" || computationSpecifier == "BY DCONST" )
    {
        mediaCorrection =std::make_shared< observation_models::ConstantReferenceCorrection >(
                        startTime, endTime, coefficients.front( ) );
    }
    else if ( computationSpecifier == "BY NRMPOW" || computationSpecifier == "BY DNRMPOW" )
    {
        mediaCorrection = std::make_shared< observation_models::PowerSeriesReferenceCorrection >(
                        startTime, endTime, coefficients );
    }
    else if ( computationSpecifier == "BY TRIG" || computationSpecifier == "BY DTRIG" )
    {
        mediaCorrection = std::make_shared< observation_models::FourierSeriesReferenceCorrection >(
                        startTime, endTime, coefficients );
    }
    else
    {
        throw std::runtime_error( "Error when creating media reference correction object: computation specifier " +
            computationSpecifier + " is not valid." );
    }

    return mediaCorrection;
}

std::vector< std::string > getGroundStationsNames( const std::string& groundStationIdentifier )
{
    std::vector< std::string > groundStations;

    if ( groundStationIdentifier == "C10" )
    {
        groundStations = { "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26" };
    }
    else if ( groundStationIdentifier == "C40" )
    {
        groundStations = { "DSS-34", "DSS-35", "DSS-36", "DSS-43", "DSS-45" };
    }
    else if ( groundStationIdentifier == "C60" )
    {
        groundStations = { "DSS-54", "DSS-55", "DSS-63", "DSS-65" };
    }
    else
    {
        groundStations = { "DSS-" + groundStationIdentifier };
    }

    return groundStations;
}

std::vector< observation_models::ObservableType > getBaseObservableTypes( const std::string& observableTypeIdentifier )
{
    std::vector< observation_models::ObservableType > baseObservableTypes;

    if ( observableTypeIdentifier == "DOPRNG" )
    {
        baseObservableTypes = { observation_models::one_way_range, observation_models::one_way_doppler };
    }
    else if ( observableTypeIdentifier == "DOPPLER" )
    {
        baseObservableTypes = { observation_models::one_way_doppler };
    }
    else if ( observableTypeIdentifier == "RANGE" )
    {
        baseObservableTypes = { observation_models::one_way_range };
    }
    else if ( observableTypeIdentifier == "VLBI" )
    {

    }
    else if ( observableTypeIdentifier == "ALL" )
    {
        baseObservableTypes = { observation_models::one_way_range, observation_models::one_way_doppler };
    }
    else
    {
        throw std::runtime_error( "Invalid CSP command observable type identifier." );
    }

    return baseObservableTypes;
}

bool compareAtmosphericCspFileStartDate( std::shared_ptr< CspRawFile > rawCspData1,
                                         std::shared_ptr< CspRawFile > rawCspData2 )
{

    std::shared_ptr< AtmosphericCorrectionCspCommand > cspCommand1 = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                rawCspData1->getCspCommands( ).at( 0 ) );
    std::shared_ptr< AtmosphericCorrectionCspCommand > cspCommand2 = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                rawCspData2->getCspCommands( ).at( 0 ) );

    if ( cspCommand1 == nullptr || cspCommand2 == nullptr )
    {
        throw std::runtime_error( "Error when comparing atmospheric CSP files: inconsistent CSP command type." );
    }

    if ( cspCommand1->startTime_ < cspCommand2->startTime_ )
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::string getSourceName(
        const std::string& sourceSpecifier,
        const int sourceId,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId,
        const std::map< int, std::string >& quasarNamePerQuasarId )
{
    std::string sourceName;

    if ( sourceSpecifier.empty( ) || sourceId == 0 )
    {
        sourceName = "";
    }
    else if ( sourceSpecifier == "SCID" )
    {
        if ( spacecraftNamePerSpacecraftId.count( sourceId ) )
        {
            sourceName = spacecraftNamePerSpacecraftId.at( sourceId );
        }
        else
        {
            // Spacecraft name selected to be "NAIF Id", which is equal to -"JPL Id" (for a spacecraft)
            // sourceName = std::to_string( - sourceId );

            throw std::runtime_error( "Error when retrieving source name from CSP atmosphere file (troposphere/ionosphere corrections):"
                                      " no spacecraft name was provided." );
        }
    }
    else if ( sourceSpecifier == "QUASAR" )
    {
        if ( quasarNamePerQuasarId.count( sourceId ) )
        {
            sourceName = quasarNamePerQuasarId.at( sourceId );
        }
        else
        {
            throw std::runtime_error( "Error when retrieving source name from CSP atmosphere file (troposphere/ionosphere corrections): "
                                      "no quasar name was provided." );
        }
    }
    else
    {
        throw std::runtime_error( "Error when retrieving source name from CSP file: invalid source specifier." );
    }

    return sourceName;
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractAtmosphericCorrection(
        std::vector< std::shared_ptr< CspRawFile > > rawCspFiles,
        const std::string& modelIdentifier,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId,
        const std::map< int, std::string >& quasarNamePerQuasarId )
{
    // Sort CSP files by start date
    std::stable_sort( rawCspFiles.begin( ), rawCspFiles.end( ), &compareAtmosphericCspFileStartDate );

    observation_models::AtmosphericCorrectionPerStationAndSpacecraftType troposphericCorrection;

    // Keys: std::pair< station, source (i.e. spacecraft or quasar) >, data type identifier
    std::map< std::pair< std::string, std::string >, std::map< std::string,
        std::shared_ptr< observation_models::TabulatedMediaReferenceCorrectionManager > > >
        correctionsPerStationsAndSourcePerType;

    for ( unsigned int i = 0; i < rawCspFiles.size( ); ++i )
    {
        for ( unsigned int j = 0; j < rawCspFiles.at( i )->getCspCommands( ).size( ); ++j )
        {
            std::shared_ptr< AtmosphericCorrectionCspCommand > cspCommand = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                    rawCspFiles.at( i )->getCspCommands( ).at( j ) );
            if ( cspCommand == nullptr )
            {
                throw std::runtime_error( "Error when creating tabulated atmospheric corrections: inconsistent CSP command type." );
            }

            if ( modelIdentifier != cspCommand->modelIdentifier_ )
            {
                continue;
            }

            std::string sourceName = getSourceName(  cspCommand->sourceSpecifier_, cspCommand->sourceId_,
                                                     spacecraftNamePerSpacecraftId, quasarNamePerQuasarId );
            if ( ( modelIdentifier == "DRY NUPART" || modelIdentifier == "WET NUPART" ) && !sourceName.empty( ) )
            {
                throw std::runtime_error( "Error when creating tabulated atmospheric corrections: invalid source name was"
                                          "created for tropospheric corrections." );
            }
            std::pair< std::string, std::string > stationsAndSource = std::make_pair(
                    cspCommand->groundStationsId_, sourceName );

            // Check whether new media correction manager object should be created
            bool createNewObject = false;
            if ( correctionsPerStationsAndSourcePerType.count( stationsAndSource ) == 0 ||
                correctionsPerStationsAndSourcePerType.at( stationsAndSource ).count( cspCommand->dataTypesIdentifier_ ) == 0 )
            {
                createNewObject = true;
            }

            if ( createNewObject )
            {
                correctionsPerStationsAndSourcePerType[ stationsAndSource ][ cspCommand->dataTypesIdentifier_ ] =
                        std::make_shared< observation_models::TabulatedMediaReferenceCorrectionManager >( );
            }

            // Create correction calculator and save it to correction manager
            correctionsPerStationsAndSourcePerType.at( stationsAndSource ).at(
                    cspCommand->dataTypesIdentifier_ )->pushReferenceCorrectionCalculator( createReferenceCorrection(
                            cspCommand->startTime_, cspCommand->endTime_,
                            cspCommand->computationCoefficients_, cspCommand->computationSpecifier_ ) );
        }
    }

    // Loop over created calculator managers and assign them to map with ground station and observation type as keys
    for ( auto stationsAndSourceIt = correctionsPerStationsAndSourcePerType.begin( );
            stationsAndSourceIt != correctionsPerStationsAndSourcePerType.end( ); ++stationsAndSourceIt )
    {
        std::vector< std::string > groundStations = getGroundStationsNames( stationsAndSourceIt->first.first );
        std::string source = stationsAndSourceIt->first.second;

        for ( auto observableIt = stationsAndSourceIt->second.begin( ); observableIt != stationsAndSourceIt->second.end( ); ++observableIt )
        {
            std::vector< observation_models::ObservableType > observableTypes = getBaseObservableTypes(
                observableIt->first );

            for ( const std::string& groundStation : groundStations )
            {
                for ( const observation_models::ObservableType& observableType : observableTypes )
                {
                    troposphericCorrection[ std::make_pair( groundStation, source ) ][ observableType ] =
                            correctionsPerStationsAndSourcePerType.at( stationsAndSourceIt->first ).at( observableIt->first );
                }
            }
        }
    }

    return troposphericCorrection;
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericDryCorrectionAdjustment(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles )
{
    std::string modelIdentifier = "DRY NUPART";
    return extractAtmosphericCorrection( rawCspFiles, modelIdentifier );
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericWetCorrectionAdjustment(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles )
{
    std::string modelIdentifier = "WET NUPART";
    return extractAtmosphericCorrection( rawCspFiles, modelIdentifier );
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractIonosphericCorrection(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId,
        const std::map< int, std::string >& quasarNamePerQuasarId )
{
    std::string modelIdentifier = "CHPART";
    return extractAtmosphericCorrection( rawCspFiles, modelIdentifier, spacecraftNamePerSpacecraftId,
                                         quasarNamePerQuasarId );
}

std::shared_ptr< CspRawFile > getDsnDefaultTroposphericSeasonalModelCspFile( )
{
    std::vector< std::shared_ptr< CspCommand > > cspCommands;

    for ( unsigned int i = 0; i < 6; ++i )
    {
        std::shared_ptr< AtmosphericCorrectionCspCommand > command = std::make_shared< AtmosphericCorrectionCspCommand >( );
        command->dataTypesIdentifier_ = "DOPRNG";
        command->computationSpecifier_ = "BY TRIG";
        command->endTime_ = TUDAT_NAN;
        command->startTime_ = command->convertTime( "72/01/01", "00:00" );

        switch ( i )
        {
            case 0:
                command->modelIdentifier_ = "WET NUPART";
                command->groundStationsId_ = "C10";
                command->computationCoefficients_ =
                        { 31557600., 0.0870, -0.0360, -0.0336, 0.0002, 0.0200, 0.0008, -0.0021, -0.0036, -0.0002 };
                break;
            case 1:
                command->modelIdentifier_ = "DRY NUPART";
                command->groundStationsId_ = "C10";
                command->computationCoefficients_ =
                        { 31557600., 2.0521, 0.0082, -0.0005, -0.0004, 0.0033, -0.0015, 0.0005, -0.0011, 0.0036 };
                break;
            case 2:
                command->modelIdentifier_ = "WET NUPART";
                command->groundStationsId_ = "C40";
                command->computationCoefficients_ =
                        { 31557600., 0.1149, 0.0255, 0.0020, 0.0010, 0.0026, 0.0036, -0.0001, 0.0007, 0.0012 };
                break;
            case 3:
                command->modelIdentifier_ = "DRY NUPART";
                command->groundStationsId_ = "C40";
                command->computationCoefficients_ =
                        { 31557600., 2.1579, -0.0032, -0.0002, 0.0012, 0.0017, -0.0043, 0.0052, 0.0016, -0.0021 };
                break;
            case 4:
                command->modelIdentifier_ = "WET NUPART";
                command->groundStationsId_ = "C60";
                command->computationCoefficients_ =
                        { 31557600., 0.1255, -0.0284, -0.0273, -0.0094, 0.0005, -0.0031, -0.0003, -0.0034, -0.0013 };
                break;
            case 5:
                command->modelIdentifier_ = "DRY NUPART";
                command->groundStationsId_ = "C60";
                command->computationCoefficients_ =
                        { 31557600., 2.1094, 0.0037, -0.0010, 0.0036, 0.0019, -0.0006, 0.0021, 0.0018, -0.0004 };
                break;
            default:
                throw std::runtime_error( "Error when getting default DSN tropospheric seasonal model." );
                break;
        }

        cspCommands.push_back( command );
    }

    return std::make_shared< CspRawFile >( cspCommands );
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericDryCorrection( )
{
    std::string modelIdentifier = "DRY NUPART";
    return extractAtmosphericCorrection(
            { getDsnDefaultTroposphericSeasonalModelCspFile( ) }, modelIdentifier );
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericWetCorrection( )
{
    std::string modelIdentifier = "WET NUPART";
    return extractAtmosphericCorrection(
            { getDsnDefaultTroposphericSeasonalModelCspFile( ) }, modelIdentifier );
}

} // namespace input_output

} // namespace tudat