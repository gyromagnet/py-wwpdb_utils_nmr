
// Generated from CnsMRLexer.g4 by ANTLR 4.13.2

#pragma once


#include "antlr4-runtime.h"




class  CnsMRLexer : public antlr4::Lexer {
public:
  enum {
    Set = 1, End = 2, Noe = 3, Analysis = 4, Assign = 5, Asymptote = 6, 
    Average = 7, Bhig = 8, Ceiling = 9, Classification = 10, CountViol = 11, 
    Cv = 12, Den = 13, Distribute = 14, Ensemble = 15, Monomers = 16, Ncount = 17, 
    Nrestraints = 18, Outd = 19, Partition = 20, Potential = 21, Predict = 22, 
    Print = 23, Raverage = 24, Threshold = 25, Reset = 26, Rswitch = 27, 
    Scale = 28, SoExponent = 29, SqConstant = 30, SqExponent = 31, SqOffset = 32, 
    Taverage = 33, Temperature = 34, Initialize = 35, Update = 36, Gamma = 37, 
    Kappa = 38, Cutoff = 39, Cuton = 40, From = 41, To = 42, Peak = 43, 
    Spectrum = 44, Volume = 45, Vol = 46, Ppm1 = 47, Ppm2 = 48, Restraints = 49, 
    Dihedral = 50, Nassign = 51, Print_any = 52, Plane = 53, Group = 54, 
    Selection = 55, Weight = 56, Harmonic = 57, Exponent = 58, Normal = 59, 
    Sanisotropy = 60, Coefficients = 61, ForceConstant = 62, Coupling = 63, 
    Carbon = 64, Expectation = 65, PhiStep = 66, PsiStep = 67, Rcoil = 68, 
    Zero = 69, Proton = 70, Observed = 71, Anisotropy = 72, Amides = 73, 
    Nitrogens = 74, Oxygens = 75, RingAtoms = 76, AlphasAndAmides = 77, 
    Error = 78, Conformation = 79, Compressed = 80, Phase = 81, Size = 82, 
    Dimensions = 83, Danisotropy = 84, OneBond = 85, AngleDb = 86, DerivFlag = 87, 
    Ncs = 88, Equivalence = 89, Sigb = 90, Flags = 91, All = 92, Around = 93, 
    Atom = 94, Attribute = 95, BondedTo = 96, ByGroup = 97, ByRes = 98, 
    Chemical = 99, Fbox = 100, Hydrogen = 101, Id = 102, Known = 103, Name = 104, 
    NONE = 105, Point = 106, Cut = 107, Previous = 108, Pseudo = 109, Residue = 110, 
    Resname = 111, Saround = 112, SegIdentifier = 113, Sfbox = 114, Store1 = 115, 
    Store2 = 116, Store3 = 117, Store4 = 118, Store5 = 119, Store6 = 120, 
    Store7 = 121, Store8 = 122, Store9 = 123, Tag = 124, Vector = 125, Do_Lp = 126, 
    Identity_Lp = 127, Show = 128, Evaluate_Lp = 129, Patch = 130, Reference = 131, 
    Nil = 132, Parameter = 133, UB = 134, Mult = 135, HBonded = 136, Improper = 137, 
    NBFix = 138, NonB = 139, VDWOff = 140, Verbose = 141, For = 142, Loop = 143, 
    Tail = 144, Head = 145, Or_op = 146, And_op = 147, Not_op = 148, Comma = 149, 
    Complex = 150, Integer = 151, Logical = 152, Real = 153, Double_quote_string = 154, 
    SHARP_COMMENT = 155, EXCLM_COMMENT = 156, SMCLN_COMMENT = 157, Simple_name = 158, 
    Simple_names = 159, Integers = 160, L_paren = 161, R_paren = 162, Colon = 163, 
    Equ_op = 164, Lt_op = 165, Gt_op = 166, Leq_op = 167, Geq_op = 168, 
    Neq_op = 169, Symbol_name = 170, SPACE = 171, ENCLOSE_COMMENT = 172, 
    SECTION_COMMENT = 173, LINE_COMMENT = 174, SET_VARIABLE = 175, Abs = 176, 
    Attr_properties = 177, Comparison_ops = 178, SPACE_AP = 179, Averaging_methods = 180, 
    Class_name_AM = 181, SPACE_AM = 182, Equ_op_CL = 183, Class_name = 184, 
    SPACE_CL = 185, Equ_op_PT = 186, Potential_types = 187, Class_name_PT = 188, 
    SPACE_PT = 189, Noe_analysis = 190, SPACE_NA = 191, Exclude = 192, Include = 193, 
    End_FL = 194, Class_name_FL = 195, Any_class = 196, SPACE_FL = 197, 
    R_paren_VE = 198, Equ_op_VE = 199, Add_op_VE = 200, Sub_op_VE = 201, 
    Mul_op_VE = 202, Div_op_VE = 203, Exp_op_VE = 204, Comma_VE = 205, Integer_VE = 206, 
    Real_VE = 207, Atom_properties_VE = 208, Abs_VE = 209, Acos_VE = 210, 
    Asin_VE = 211, Cos_VE = 212, Decode_VE = 213, Encode_VE = 214, Exp_VE = 215, 
    Gauss_VE = 216, Heavy_VE = 217, Int_VE = 218, Log10_VE = 219, Log_VE = 220, 
    Max_VE = 221, Maxw_VE = 222, Min_VE = 223, Mod_VE = 224, Norm_VE = 225, 
    Random_VE = 226, Sign_VE = 227, Sin_VE = 228, Sqrt_VE = 229, Tan_VE = 230, 
    Symbol_name_VE = 231, Simple_name_VE = 232, Double_quote_string_VE = 233, 
    SPACE_VE = 234, L_paren_VF = 235, SPACE_VF = 236, L_paren_VS = 237, 
    R_paren_VS = 238, Average_VS = 239, Element_VS = 240, Max_VS = 241, 
    Min_VS = 242, Norm_VS = 243, Rms_VS = 244, Sum_VS = 245, Atom_properties_VS = 246, 
    SPACE_VS = 247, L_paren_CF = 248, R_paren_CF = 249, In_CF = 250, Integer_CF = 251, 
    Real_CF = 252, Symbol_name_CF = 253, Simple_name_CF = 254, SPACE_CF = 255, 
    COMMENT_CF = 256, Simple_name_LL = 257, SPACE_LL = 258
  };

  enum {
    ATTR_MODE = 1, AVER_MODE = 2, CLAS_MODE = 3, POTE_MODE = 4, ANAL_MODE = 5, 
    FLAG_MODE = 6, VECTOR_EXPR_MODE = 7, VECTOR_FUNC_MODE = 8, VECTOR_SHOW_MODE = 9, 
    CTL_FOR_MODE = 10, LOOP_LABEL_MODE = 11
  };

  explicit CnsMRLexer(antlr4::CharStream *input);

  ~CnsMRLexer() override;


  std::string getGrammarFileName() const override;

  const std::vector<std::string>& getRuleNames() const override;

  const std::vector<std::string>& getChannelNames() const override;

  const std::vector<std::string>& getModeNames() const override;

  const antlr4::dfa::Vocabulary& getVocabulary() const override;

  antlr4::atn::SerializedATNView getSerializedATN() const override;

  const antlr4::atn::ATN& getATN() const override;

  // By default the static state used to implement the lexer is lazily initialized during the first
  // call to the constructor. You can call this function if you wish to initialize the static state
  // ahead of time.
  static void initialize();

private:

  // Individual action functions triggered by action() above.

  // Individual semantic predicate functions triggered by sempred() above.

};

