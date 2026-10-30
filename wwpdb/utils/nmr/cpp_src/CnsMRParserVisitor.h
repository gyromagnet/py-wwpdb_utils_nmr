
// Generated from CnsMRParser.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"
#include "CnsMRParser.h"



/**
 * This class defines an abstract visitor for a parse tree
 * produced by CnsMRParser.
 */
class  CnsMRParserVisitor : public antlr4::tree::AbstractParseTreeVisitor {
public:

  /**
   * Visit parse trees produced by CnsMRParser.
   */
    virtual std::any visitCns_mr(CnsMRParser::Cns_mrContext *context) = 0;

    virtual std::any visitDistance_restraint(CnsMRParser::Distance_restraintContext *context) = 0;

    virtual std::any visitDihedral_angle_restraint(CnsMRParser::Dihedral_angle_restraintContext *context) = 0;

    virtual std::any visitPlane_restraint(CnsMRParser::Plane_restraintContext *context) = 0;

    virtual std::any visitHarmonic_restraint(CnsMRParser::Harmonic_restraintContext *context) = 0;

    virtual std::any visitRdc_restraint(CnsMRParser::Rdc_restraintContext *context) = 0;

    virtual std::any visitCoupling_restraint(CnsMRParser::Coupling_restraintContext *context) = 0;

    virtual std::any visitCarbon_shift_restraint(CnsMRParser::Carbon_shift_restraintContext *context) = 0;

    virtual std::any visitProton_shift_restraint(CnsMRParser::Proton_shift_restraintContext *context) = 0;

    virtual std::any visitConformation_db_restraint(CnsMRParser::Conformation_db_restraintContext *context) = 0;

    virtual std::any visitDiffusion_anisotropy_restraint(CnsMRParser::Diffusion_anisotropy_restraintContext *context) = 0;

    virtual std::any visitOne_bond_coupling_restraint(CnsMRParser::One_bond_coupling_restraintContext *context) = 0;

    virtual std::any visitAngle_db_restraint(CnsMRParser::Angle_db_restraintContext *context) = 0;

    virtual std::any visitNoe_statement(CnsMRParser::Noe_statementContext *context) = 0;

    virtual std::any visitNoe_assign(CnsMRParser::Noe_assignContext *context) = 0;

    virtual std::any visitPredict_statement(CnsMRParser::Predict_statementContext *context) = 0;

    virtual std::any visitNoe_annotation(CnsMRParser::Noe_annotationContext *context) = 0;

    virtual std::any visitDihedral_statement(CnsMRParser::Dihedral_statementContext *context) = 0;

    virtual std::any visitDihedral_assign(CnsMRParser::Dihedral_assignContext *context) = 0;

    virtual std::any visitPlane_statement(CnsMRParser::Plane_statementContext *context) = 0;

    virtual std::any visitPlane_group(CnsMRParser::Plane_groupContext *context) = 0;

    virtual std::any visitGroup_statement(CnsMRParser::Group_statementContext *context) = 0;

    virtual std::any visitHarmonic_statement(CnsMRParser::Harmonic_statementContext *context) = 0;

    virtual std::any visitHarmonic_assign(CnsMRParser::Harmonic_assignContext *context) = 0;

    virtual std::any visitSani_statement(CnsMRParser::Sani_statementContext *context) = 0;

    virtual std::any visitSani_assign(CnsMRParser::Sani_assignContext *context) = 0;

    virtual std::any visitCoupling_statement(CnsMRParser::Coupling_statementContext *context) = 0;

    virtual std::any visitCoup_assign(CnsMRParser::Coup_assignContext *context) = 0;

    virtual std::any visitCarbon_shift_statement(CnsMRParser::Carbon_shift_statementContext *context) = 0;

    virtual std::any visitCarbon_shift_assign(CnsMRParser::Carbon_shift_assignContext *context) = 0;

    virtual std::any visitCarbon_shift_rcoil(CnsMRParser::Carbon_shift_rcoilContext *context) = 0;

    virtual std::any visitProton_shift_statement(CnsMRParser::Proton_shift_statementContext *context) = 0;

    virtual std::any visitObserved(CnsMRParser::ObservedContext *context) = 0;

    virtual std::any visitProton_shift_rcoil(CnsMRParser::Proton_shift_rcoilContext *context) = 0;

    virtual std::any visitProton_shift_anisotropy(CnsMRParser::Proton_shift_anisotropyContext *context) = 0;

    virtual std::any visitProton_shift_amides(CnsMRParser::Proton_shift_amidesContext *context) = 0;

    virtual std::any visitProton_shift_carbons(CnsMRParser::Proton_shift_carbonsContext *context) = 0;

    virtual std::any visitProton_shift_nitrogens(CnsMRParser::Proton_shift_nitrogensContext *context) = 0;

    virtual std::any visitProton_shift_oxygens(CnsMRParser::Proton_shift_oxygensContext *context) = 0;

    virtual std::any visitProton_shift_ring_atoms(CnsMRParser::Proton_shift_ring_atomsContext *context) = 0;

    virtual std::any visitProton_shift_alphas_and_amides(CnsMRParser::Proton_shift_alphas_and_amidesContext *context) = 0;

    virtual std::any visitConformation_statement(CnsMRParser::Conformation_statementContext *context) = 0;

    virtual std::any visitConf_assign(CnsMRParser::Conf_assignContext *context) = 0;

    virtual std::any visitDiffusion_statement(CnsMRParser::Diffusion_statementContext *context) = 0;

    virtual std::any visitDani_assign(CnsMRParser::Dani_assignContext *context) = 0;

    virtual std::any visitOne_bond_coupling_statement(CnsMRParser::One_bond_coupling_statementContext *context) = 0;

    virtual std::any visitOne_bond_assign(CnsMRParser::One_bond_assignContext *context) = 0;

    virtual std::any visitAngle_db_statement(CnsMRParser::Angle_db_statementContext *context) = 0;

    virtual std::any visitAngle_db_assign(CnsMRParser::Angle_db_assignContext *context) = 0;

    virtual std::any visitNcs_restraint(CnsMRParser::Ncs_restraintContext *context) = 0;

    virtual std::any visitNcs_statement(CnsMRParser::Ncs_statementContext *context) = 0;

    virtual std::any visitNcs_group_statement(CnsMRParser::Ncs_group_statementContext *context) = 0;

    virtual std::any visitSelection(CnsMRParser::SelectionContext *context) = 0;

    virtual std::any visitSelection_expression(CnsMRParser::Selection_expressionContext *context) = 0;

    virtual std::any visitTerm(CnsMRParser::TermContext *context) = 0;

    virtual std::any visitFactor(CnsMRParser::FactorContext *context) = 0;

    virtual std::any visitNumber(CnsMRParser::NumberContext *context) = 0;

    virtual std::any visitNumber_f(CnsMRParser::Number_fContext *context) = 0;

    virtual std::any visitNumber_s(CnsMRParser::Number_sContext *context) = 0;

    virtual std::any visitNumber_a(CnsMRParser::Number_aContext *context) = 0;

    virtual std::any visitClassification(CnsMRParser::ClassificationContext *context) = 0;

    virtual std::any visitClass_name(CnsMRParser::Class_nameContext *context) = 0;

    virtual std::any visitFlag_statement(CnsMRParser::Flag_statementContext *context) = 0;

    virtual std::any visitVector_statement(CnsMRParser::Vector_statementContext *context) = 0;

    virtual std::any visitVector_mode(CnsMRParser::Vector_modeContext *context) = 0;

    virtual std::any visitVector_expression(CnsMRParser::Vector_expressionContext *context) = 0;

    virtual std::any visitVector_operation(CnsMRParser::Vector_operationContext *context) = 0;

    virtual std::any visitVflc(CnsMRParser::VflcContext *context) = 0;

    virtual std::any visitVector_func_call(CnsMRParser::Vector_func_callContext *context) = 0;

    virtual std::any visitVector_show_property(CnsMRParser::Vector_show_propertyContext *context) = 0;

    virtual std::any visitEvaluate_statement(CnsMRParser::Evaluate_statementContext *context) = 0;

    virtual std::any visitEvaluate_operation(CnsMRParser::Evaluate_operationContext *context) = 0;

    virtual std::any visitPatch_statement(CnsMRParser::Patch_statementContext *context) = 0;

    virtual std::any visitParameter_setting(CnsMRParser::Parameter_settingContext *context) = 0;

    virtual std::any visitParameter_statement(CnsMRParser::Parameter_statementContext *context) = 0;

    virtual std::any visitNoe_assign_loop(CnsMRParser::Noe_assign_loopContext *context) = 0;

    virtual std::any visitDihedral_assign_loop(CnsMRParser::Dihedral_assign_loopContext *context) = 0;

    virtual std::any visitSani_assign_loop(CnsMRParser::Sani_assign_loopContext *context) = 0;

    virtual std::any visitCoup_assign_loop(CnsMRParser::Coup_assign_loopContext *context) = 0;

    virtual std::any visitCarbon_shift_assign_loop(CnsMRParser::Carbon_shift_assign_loopContext *context) = 0;

    virtual std::any visitPlane_group_loop(CnsMRParser::Plane_group_loopContext *context) = 0;


};

