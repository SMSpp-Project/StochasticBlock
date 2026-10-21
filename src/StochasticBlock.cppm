/*--------------------------------------------------------------------------*/
/*----------------------- File StochasticBlock.cppm ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.stochasticblock, which exports every
 * name that the headers of StochasticBlock declare at namespace scope,
 * together with the modules it depends on; a translation unit can then write
 *
 *     import smspp.stochasticblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "DiscreteScenarioSet.h"
#include "IndependentMultiStageScenarioGenerator.h"
#include "MultiStageDiscreteScenarioSet.h"
#include "ScenarioGenerator.h"
#include "ScenarioReductionBlock.h"
#include "StochasticBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.stochasticblock;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // DiscreteScenarioSet.h
 using SMSpp_di_unipi_it::DiscreteScenarioSet;

 // IndependentMultiStageScenarioGenerator.h
 using SMSpp_di_unipi_it::IndependentMultiStageScenarioGenerator;

 // MultiStageDiscreteScenarioSet.h
 using SMSpp_di_unipi_it::MultiStageDiscreteScenarioSet;

 // ScenarioGenerator.h
 using SMSpp_di_unipi_it::MultiStageScenarioGenerator;
 using SMSpp_di_unipi_it::ScenarioGenerator;

 // ScenarioReductionBlock.h
 using SMSpp_di_unipi_it::ScenarioReductionBlock;
 using SMSpp_di_unipi_it::ScenarioReductionBlockSolution;

 // StochasticBlock.h
 using SMSpp_di_unipi_it::StochasticBlock;
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File StochasticBlock.cppm ----------------------*/
/*--------------------------------------------------------------------------*/
