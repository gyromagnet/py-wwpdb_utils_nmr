/*
 CNS MR (Magnetic Restraint) lexer grammar for ANTLR v4.10 or later
 Copyright 2022 Masashi Yokochi

you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

lexer grammar CnsMRLexer;

options { caseInsensitive=true; }

Set:			'SET';
End:			'END';

/* CNS: Distance restraints - Syntax - noe
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Noe:			'NOE';					// Noe { noe_statement } End

// noe_statement
Analysis:		'ANAL' 'Y'? 'S'? 'I'? 'S'? -> pushMode(ANAL_MODE);		// = Noe_analysis
Assign:			'ASSI' 'G'? 'N'?;			// selection selection Real Real Real [ Or_op selection selection ... ]
Asymptote:		'ASYM' 'P'? 'T'? 'O'? 'T'? 'E'?;	// Class_names Real
Average:		'AVER' 'A'? 'G'? (('I'? 'N'? 'G'?) | 'E'?) -> pushMode(AVER_MODE);	// Class_names Noe_avr_methods
Bhig:			'BHIG';					// Class_names Real
Ceiling:		'CEIL' 'I'? 'N'? 'G'?;			// = Real
Classification:		'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'? -> pushMode(CLAS_MODE);	// Class_name
CountViol:		'COUN' 'T'? 'V'? 'I'? 'O'? 'L'?;	// Class_name
Cv:			'CV';					// = Integer
Den:			'DEN';					// Initialize | ( Update Gamma = Real Kappa = Real )
Distribute:		'DIST' 'R'? 'I'? 'B'? 'U'? 'T'? 'E'?;	// Class_name Class_name Real
Ensemble:		'ENSE' 'M'? 'B'? 'L'? 'E'?;		// Ensemble { ensemble_statement } End
Monomers:		'MONO' 'M'? 'E'? 'R'? 'S'?;		// Class_names Integer
Ncount:			'NCOU' 'N'? 'T'?;			// Class_names Integer
Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
Outd:			'OUTD';
Partition:		'PART' 'I'? 'T'? 'I'? 'O'? 'N'?;	// = Integer
Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// Class_names Noe_potential
Predict:		'PRED' 'I'? 'C'? 'T'?;			// { predict_statement } End
Print:			'PRIN' 'T'?;
Raverage:		'RAVE' 'R'? 'A'? 'G'? 'E'?;		// Class_name Raverage_statement End
Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// = Real
Reset:			'RESE' 'T'?;
Rswitch:		'RSWI' 'T'? 'C'? 'H'?;			// Class_names Real
Scale:			'SCAL' 'E'?;				// Class_names Real
SoExponent:		'SOEX' 'P'? 'O'? 'N'? 'E'? 'N'? 'T'?;	// Class_names Real
SqConstant:		'SQCO' 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;			// Class_names Real
SqExponent:		'SQEX' 'P'? 'O'? 'N'? 'E'? 'N'? 'T'?;	// Class_names Real
SqOffset:		'SQOF' 'F'? 'S'? 'E'? 'T'?;		// Class_names Real
Taverage:		'TAVE' 'R'? 'A'? 'G'? 'E'?;		// Class_name Taverage_statement End
Temperature:		'TEMP' 'E'? 'R'? 'A'? 'T'? 'U'? 'R'? 'E'?;			// = Real

// NOE analysis
//Noe_analysis:		'CURR' 'E'? 'N'? T? | 'TAVE' 'R'? 'A'? 'G'? 'E'? | 'RAVE' 'R'? 'A'? 'G'? 'E'?;

Initialize:		'INIT' 'I'? 'A'? 'L'? 'I'? 'Z'? 'E'?;
Update:			'UPDA' 'T'? 'E'?;			// Gamma = Real Kappa = Real
Gamma:			'GAMM' 'A'?;
Kappa:			'KAPP' 'A'?;

// NOE averaging methods
//Noe_avr_methods:	'R-6' | 'R-3' | 'SUM' | 'CENT' 'E'? 'R'?;

// NOE potential statement
//Noe_potential:	'BIHA' 'R'? 'M'? 'O'? 'N'? 'I'? 'C'? | 'LOGN' 'O'? 'R'? 'M'? 'A'? 'L'? | 'SQUA' 'R'? 'E'? '-'? 'W'? 'E'? 'L'? 'L'? | 'SOFT' '-'? 'S'? 'Q'? 'U'? 'A'? 'R'? 'E'? | 'SYMM' 'E'? 'T'? 'R'? 'Y'? | 'HIGH' | '3DPO';

// Predict statement
Cutoff:			'CUTOFF';				// = Real
Cuton:			'CUTON';				// = Real
From:			'FROM';					// = selection
To:			'TO';					// = selection

// 3rd party software extensions for NOE assign clause
Peak:			'PEAK';					// = Integer
Spectrum:		'SPECTRUM';				// = Integer
//Weight:		'WEIGHT';				// = Real
Volume:			'VOLUME';				// = Real
Vol:			'VOL';					// = Real
Ppm1:			'PPM1';					// = Real
Ppm2:			'PPM2';					// = Real
//Cv:			'CV';					// = Integer

/* CNS: Dihedral angle restraints - Syntax - restranits/dihedral
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Restraints:		'REST' 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;	// Dihedral
Dihedral:		'DIHE' 'D'? 'R'? 'A'? 'L'?;		// Dihedral { dihedral_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection Real Real Real Integer
//Cv:			'CV';					// = Integer
Nassign:		'NASS' 'I'? 'G'? 'N'?;			// = Integer
//Partition:		'PART' 'I'? 'T'? 'I'? 'O'? 'N'?;	// = Integer
//Reset:		'RESE' 'T'?;
//Scale:		'SCAL' 'E'?;				// = Real
Print_any:		'?';

/* CNS: Plane restraints - Syntax - restraints/plane
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
//Restraints:		'REST' 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;	// Plane
Plane:			'PLAN' ('E' | 'A')? 'R'?;		// Plane { planar_statement } End

// planar_statement
Group:			'GROU' 'P'?;				// Group { planar_group_statement } End
//Initialize:		'INIT' 'I'? 'A'? 'L'? 'I'? 'Z'? 'E'?;
//Print_any:		'?';

// group_statement
Selection:		'SELE' 'C'? 'T'? 'I'? 'O'? 'N'? '='?;	// = selection
Weight:			'WEIG' 'H'? 'T'?;			// = Real

/* CNS: Plane restraints - Syntax - restraints/harmonic
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
//Restraints:		'REST' 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;	// Harmonic
Harmonic:		'HARM' 'O'? 'N'? 'I'? 'C'?;		// Harmonic { harmonic_statement } End

// harmonic_stetement
Exponent:		'EXPO' 'N'? 'E'? 'N'? 'T'?;		// = Integer
Normal:			'NORM' 'A'? 'L'?;			// = vector_3d

/* CNS: Suscetibility anisotropy restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Sanisotropy:		'SANI' 'S'? 'O'? 'T'? 'R'? 'O'? 'P'? 'Y'?;			// Sanisotropy { sani_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection selection selection Real Real
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
Coefficients:		'COEF' 'F'? 'I'? 'C'? 'I'? 'E'? 'N'? 'T'? 'S'?;			// Real Real Real
ForceConstant:		'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// = Real
//Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// = Rdc_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real
//Reset:		'RESE' 'T'?;

// RDC potential statement
//Rdc_potential:	'SQUA' 'R'? 'E'? '-'? 'W'? 'E'? 'L'? 'L'? | 'HARM' 'O'? 'N'? 'I'? 'C'?;

/* CNS: Scalar J-coupling restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Coupling:		'COUP' 'L'? 'I'? 'N'? 'G'?;		// Coupling { coupling_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection [ selection selection selection selection ] Real Real [ Real Real ]
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
//Coefficients:		'COEF' 'F'? 'I'? 'C'? 'I'? 'E'? 'N'? 'T'? 'S'?;			// Real Real Real Real
//Cv:			'CV';					// = Integer
//ForceConstant:	'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// Real [ Real ]
//Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
//Partition:		'PART' 'I'? 'T'? 'I'? 'O'? 'N'?;	// = Integer
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// = Coupling_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real ( All | Class = Class_name )
//Reset:		'RESE' 'T'?;

//Coupling_potential:	Rdc_potential | 'MULT' 'I'? 'P'? 'L'? 'E'?;

/* CNS: Carbon chemical shift restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Carbon:			'CARB' 'O'? 'N'?;			// Carbon { carbon_shift_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection selection Real Real
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
Expectation:		'EXPE' 'C'? 'T'? 'A'? 'T'? 'I'? 'O'? 'N'?;			// Integer Integer Real Real Real Real
//ForceConstant:	'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// = Real
//Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
PhiStep:		'PHIS' 'T'? 'E'? 'P'?;			// = Real
PsiStep:		'PSIS' 'T'? 'E'? 'P'?;			// = Real
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// = Rdc_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real
Rcoil:			'RCOI' 'L'?;				// selection Real Real
//Reset:		'RESE' 'T'?;
Zero:			'ZERO';

/* CNS: Proton chemical shift restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Proton:			'PROTON' 'S'? 'H'? 'I'? 'F'? 'T'? 'S'?;	// Proton { proton_shift_statement } End
Observed:		'OBSE' 'R'? 'V'? 'E'? 'D'?;		// selection [ selection ] Real [ Real ]
//Rcoil:		'RCOI' 'L'?;				// selection Real
Anisotropy:		'ANIS' 'O'? 'T'? 'R'? 'O'? 'P'? 'Y'?;	// selection selection selection Co_or_Cn Logical? SC_or_BB
Amides:			'AMID' 'E'? 'S'?;			// selection
//Carbons:		'CARB' 'O'? 'N'? 'S'?;			// selection
Nitrogens:		'NITR' 'O'? 'G'? 'E'? 'N'? 'S'?;	// selection
Oxygens:		'OXYG' 'E'? 'N'? 'S'?;			// selection
RingAtoms:		'RING' 'A'? 'T'? 'O'? 'M'? 'S'?;	// Ring_resname selection selection selection selection selection [ selection ]
AlphasAndAmides:	'ALPH' 'A'? 'S'? 'A'? 'N'? 'D'? 'A'? 'M'? 'I'? 'D'? 'E'? 'S'?;	// selection
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
Error:			'ERRO' 'R'?;				// Real
//ForceConstant:	'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// Real [ Real ]
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// Coupling_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real ( All | Class = Class_name ) Rmsd_or_Not
//Reset:		'RESE' 'T'?;

//CO_or_CN:		'CO' | 'CN';
//SC_or_BB:		'SC' | 'BB';
//Ring_resname:		'PHE' | 'TYR' | 'HIS' | 'TRP' ('5' | '6') | 'ADE' ('5' | '6') | 'GUA' ('5' | '6') | 'THY' | 'CYT' | 'URA';
//Rmsd_or_Not:		'RMSD' | 'NORM' 'S'? 'D'?;

/* CNS: Conformation database restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Conformation:		'CONF' 'O'? 'R'? 'M'? 'A'? 'T'? 'I'? 'O'? 'N'?;			// Conformation { conformation_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection [ selection selection selection selection ] [ selection selection selection selection ] [ selection selection selection selection ]
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
Compressed:		'COMP' 'R'? 'E'? 'S'? 'S'? 'E'? 'D'?;
//Expectation:		'EXPE' 'C'? 'T'? 'A'? 'T'? 'I'? 'O'? 'N'?;			// Integer [ Integer ] [ Integer ] [ Integer ] Real
//Error:		'ERRO' 'R'?;				// = Real
//ForceConstant:	'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// = Real
//Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
Phase:			'PHAS' 'E'?;				// Integer Integer Integer [ Integer Integer Integer ] (Integer Integer Integer ] (Integer Integer Integer ]
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// = Rdc_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real ( All | Class = Class_name )
//Reset:		'RESE' 'T'?;
Size:			'SIZE';					// Dimensions Integer [ Integer ] [ Integer ] [ Integer ]
//Zero:			'ZERO';

Dimensions:		'ONED' | 'TWOD' | 'THREED' | 'FOURD';

/* CNS: Diffusion anisotropy restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Danisotropy:		'DANI' 'S'? 'O'? 'T'? 'R'? 'O'? 'P'? 'Y'?;			// Danisotropy { diffusion_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection selection selection Real Real
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
//Coefficients:		'COEF' 'F'? 'I'? 'C'? 'I'? 'E'? 'N'? 'T'? 'S'?;			// Real Real Real Real Real
//ForceConstant:	'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// = Real
//Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// = Rdc_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real
//Reset:		'RESE' 'T'?;

/* CNS: One-bond coupling restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
OneBond:		'ONEB' 'O'? 'N'? 'D'?;			// OneBond { one_bond_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection selection selection selection selection Real Real
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
//Coefficients:		'COEF' 'F'? 'I'? 'C'? 'I'? 'E'? 'N'? 'T'? 'S'?;			// Real Real Real Real Real Real Real
//ForceConstant:	'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// = Real
//Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// = Rdc_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real
//Reset:		'RESE' 'T'?;

/* CNS: Angle database restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
AngleDb:		'ANGL' 'E'? 'D'? 'A'? 'T'? 'A'? 'B'? 'A'? 'S'? 'E'?;		// AngleDb { bond_angle_statement } End
//Assign:		'ASSI' 'G'? 'N'?;			// selection selection selection selection selection selection selection selection selection selection selection [ selection ]
//Classification:	'CLAS' 'S'? 'I'? 'F'? 'I'? 'C'? 'A'? 'T'? 'I'? 'O'? 'N'?;	// = Class_name
DerivFlag:		'DERI' 'V'? 'F'? 'L'? 'A'? 'G'?;	// On_or_Off
//Expectation:		'EXPE' 'C'? 'T'? 'A'? 'T'? 'I'? 'O'? 'N'?;			// Integer Integer Real
//Error:		'ERRO' 'R'?;				// = Real
//ForceConstant:	'FORC' 'E'? 'C'? 'O'? 'N'? 'S'? 'T'? 'A'? 'N'? 'T'?;		// = Real
//Nrestraints:		'NRES' 'T'? 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;			// = Integer
//Potential:		'POTE' 'N'? 'T'? 'I'? 'A'? 'L'? -> pushMode(POTE_MODE);		// = Rdc_potential
//Print:		'PRIN' 'T'?;
//Threshold:		'THRE' 'S'? 'H'? 'O'? 'L'? 'D'?;	// Real ( All | Class = Class_name )
//Reset:		'RESE' 'T'?;
//Size:			'SIZE';					// Angle_or_Dihedral Integer Integer
//Zero:			'ZERO';

//On_or_Off:		'ON' | 'OFF';
//Angle_or_Dihedral:	'ANGL' 'E'? | 'DIHE' 'D'? 'R'? 'A'? 'L'?;

/* CNS: NCS restraints - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Ncs:			'NCS';
//Restraints:		'REST' 'R'? 'A'? 'I'? 'N'? 'T'? 'S'?;	// Ncs

//Group:		'GROU' 'P'?;				// Group { ncs_group_statement } End
//Initialize:		'INIT' 'I'? 'A'? 'L'? 'I'? 'Z'? 'E'?;
//Print_any:		'?';

Equivalence:		'EQUI' 'V'? 'A'? 'L'? 'E'? 'N'? 'C'? 'E'? '='?;
								// = selection
Sigb:			'SIGB';					// = Real
//Weight:		'WEIG' 'H'? 'T'?;			// = Real

/* CNS: Flags - Syntax
 See alos https://nmr.cit.nih.gov/xplor-nih/xplorMan/node125.html (compatible with XPLOR-NIH)
*/
Flags:			'FLAG' 'S'? -> pushMode(FLAG_MODE);	// Flags { flag_statement } End

/* Atom selection - Syntax - identity/atom-selection
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
All:			'ALL';
Around:			'AROU' 'N'? 'D'?;			// Real (factor as subject)
Atom:			'ATOM';					// Segment_names Residue_numbers Atom_names
Attribute:		'ATTR' 'I'? 'B'? 'U'? 'T'? 'E'? -> pushMode(ATTR_MODE);		// Abs? Attr_property Comparison_ops Real
BondedTo:		'BOND' 'E'? 'D'? 'T'? 'O'?;		// factor
ByGroup:		'BYGR' 'O'? 'U'? 'P'?;			// factor
ByRes:			'BYRE' 'S'?;				// factor
Chemical:		'CHEM' 'I'? 'C'? 'A'? 'L'?;		// Atom_types | Atom_type [ : Atom_type ]
Fbox:			'FBOX';					// real real real real real real
Hydrogen:		'HYDR' 'O'? 'G'? 'E'? 'N'?;
Id:			'ID';					// Integer
Known:			'KNOW' 'N'?;
Name:			'NAME';					// Atom_names | Atom_name [ : Atom_name ]
NONE:			'NONE';
//Not_op:		'NOT';					// factor
Point:			'POIN' 'T'?;				// vector_3d cut Real
Cut:			'CUT';
Previous:		'PREV' 'I'? 'O'? 'U'? 'S'?;
Pseudo:			'PSEU' 'D'? 'O'?;
Residue:		'RESI' 'D'? 'U'? 'E'?;			// Residue_numbers | Residue_number [ : Residue_number ]
Resname:		'RESN' 'A'? 'M'? 'E'?;			// Residue_names | Residue_name [ : Residue_name ]
Saround:		'SARO' 'U'? 'N'? 'D'?;			// Real (factor as subject)
SegIdentifier:		'SEGI' 'D'? 'E'? 'N'? 'T'? 'I'? 'F'? 'I'? 'E'? 'R'?;		// Segment_names | Segment_name [ : Segment_name ]
Sfbox:			'SFBOX';				// real real real real real real
Store1:			'STORE1';
Store2:			'STORE2';
Store3:			'STORE3';
Store4:			'STORE4';
Store5:			'STORE5';
Store6:			'STORE6';
Store7:			'STORE7';
Store8:			'STORE8';
Store9:			'STORE9';
Tag:			'TAG';

/* Vector statement - Syntax
 See also https://nmr.cit.nih.gov/xplor-nih/xplorMan/node42.html (compatible with XPLOR-NIH)
*/
Vector:			'VECT' 'O'? 'R'?;			// vector_mode vector_expression selection

Do_Lp:			'DO' ' '* L_paren -> pushMode(VECTOR_EXPR_MODE);
Identity_Lp:		'IDEN' 'T'? 'I'? 'T'? 'Y'? ' '* L_paren -> pushMode(VECTOR_EXPR_MODE);
Show:			'SHOW' -> pushMode(VECTOR_SHOW_MODE);	// Vector_show_property

/* CNS: Global statement/Evaluate statement - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Evaluate_Lp:		'EVAL' 'U'? 'A'? 'T'? 'E'? ' '* L_paren -> pushMode(VECTOR_EXPR_MODE);
								// ( evaluate_statement )

/* CNS: Patch molecular topology database - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Patch:			'PATC' 'H'?;				// ( patch_statement )
Reference:		'REFE' 'R'? 'E'? 'N'? 'C'? 'E'?;
Nil:			'NIL';

/* CNS: Parameter statement - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Parameter:		'PARA' 'M'? 'E'? 'T'? 'E'? 'R'?;	// ( parameter_statement )
UB:			'UB';
Mult:			'MULT';
HBonded:		'HBON' 'D'? 'E'? 'D'?;
Improper:		'IMPR' 'O'? 'P'? 'E'? 'R'?;
NBFix:			'NBFI' 'X'?;
NonB:			'NONB';
VDWOff:			'VDWO' 'F'? 'F'?;
Verbose:		'VERB' 'O'? 'S'? 'E'?;

/* Control statement - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
For:			'FOR' -> pushMode(CTL_FOR_MODE);	// Symbol_name In ( Words ) Loop Loop_label { statements } End Loop Loop_label
Loop:			'LOOP' -> pushMode(LOOP_LABEL_MODE);

/* Three-dimentional vectors - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Tail:			'TAIL' '='?;
Head:			'HEAD' '='?;

// Logical operations
Or_op:			'OR';
And_op:			'AND';
Not_op:			'NOT';

/* Numbers and Strings - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
Comma:			',';
Complex:		L_paren Real Comma Real R_paren;
Integer:		'-'? DECIMAL;
Logical:		'TRUE' | 'FALSE' | 'ON' | 'OFF';
Real:			('+' | '-')? (DECIMAL | DEC_DOT_DEC) ('E' ('+' | '-')? DECIMAL)?;
Double_quote_string:	'"' ~["\r\n]* '"';
fragment DEC_DOT_DEC:	(DECIMAL '.' DECIMAL?) | ('.' DECIMAL);
fragment DEC_DIGIT:	[0-9];
fragment DECIMAL:	DEC_DIGIT+;

SHARP_COMMENT:		'#'+ ~[\r\n]* '#'* ~[\r\n]* -> channel(HIDDEN);
EXCLM_COMMENT:		'!'+ ~[\r\n]* '!'* ~[\r\n]* -> channel(HIDDEN);
SMCLN_COMMENT:		';'+ ~[\r\n]* ';'* ~[\r\n]* -> channel(HIDDEN);

Simple_name:		SIMPLE_NAME;
Simple_names:		(WILDCARD | WILDCARD* SIMPLE_NAME WILDCARD+ (SIMPLE_NAME | WILDCARD)*) POST_WC_CHAR*;
Integers:		(WILDCARD | WILDCARD* Integer WILDCARD+) DEC_DIGIT*;

/* Wildcard - Syntax
 See also https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/syntax_manual/frame.html
*/
fragment WILDCARD:	'*' | '%' | '#' | '+';

fragment ALPHA:		[A-Z];
fragment ALPHA_NUM:	ALPHA | DEC_DIGIT;
fragment START_CHAR:	ALPHA_NUM | '_';
fragment NAME_CHAR:	START_CHAR | '\'' | '-' | '+' | '.' | '"' | '/';
//fragment ATM_NAME_CHAR:	ALPHA_NUM | '\'';
//fragment ATM_TYPE_CHAR:	ALPHA_NUM | '-' | '+';
fragment SIMPLE_NAME:	START_CHAR NAME_CHAR*;
fragment POST_WC_CHAR:	DEC_DIGIT | '\'' | 'P';
fragment SYMBOL_NAME:	'$' START_CHAR+;

L_paren:		'(';
R_paren:		')';
Colon:			':';
Equ_op:			'=';
Lt_op:			'<';
Gt_op:			'>';
Leq_op:			'<=';
Geq_op:			'>=';
Neq_op:			'#';

Symbol_name:		SYMBOL_NAME;

SPACE:			[ \t\r\n]+ -> skip;
ENCLOSE_COMMENT:	'{' (ENCLOSE_COMMENT | .)*? '}' -> channel(HIDDEN);
SECTION_COMMENT:	('#' | '!' | ';' | '\\' | '/' '/'+ | '*' '*'+ | '-' '-'+ | '+' '+'+ | '=' '='+ | '>' '>'+ | 'REMARK' 'S'?) ' '* [\r\n]+ -> channel(HIDDEN);
LINE_COMMENT:		('#' | '!' | ';' | '\\' | '/' '/'+ | '*' '*'+ | '-' '-'+ | '+' '+'+ | '=' '='+ | '>' '>'+ | 'REMARK' 'S'?) ~[\r\n]* -> channel(HIDDEN);
SET_VARIABLE:		Set ([\r\n]*)? ~[\r\n]* ([\r\n]*)? End -> channel(HIDDEN);

mode ATTR_MODE; // Inside of Attribute tag

// Attribute properties
Abs:			'ABS';
Attr_properties:	('B' | 'BCOM' 'P'? | 'CHAR' 'G'? 'E'? | 'DX' | 'DY' | 'DZ' | 'FBET' 'A'? | 'HARM' 'O'? 'N'? 'I'? 'C'? 'S'? | 'MASS' | 'Q' | 'QCOM' 'P'? | 'REFX' | 'REFY' | 'REFZ' | 'RMSD' | 'VX' | 'VY' | 'VZ' | 'X' | 'XCOM' 'P'? | 'Y' | 'YCOM' 'P'? | 'Z' | 'ZCOM' 'P'? | 'STORE1' | 'STORE2' | 'STORE3' | 'STORE4' | 'STORE5' | 'STORE6' | 'STORE7' | 'STORE8' | 'STORE9' | 'SCATTER_A1' | 'SCATTER_A2' | 'SCATTER_A3' | 'SCATTER_A4' | 'SCATTER_B1' | 'SCATTER_B2' | 'SCATTER_B3' | 'SCATTER_B4' | 'SCATTER_C' | 'SCATTER_FP' | 'SCATTER_FDP');
Comparison_ops:		(Equ_op | Lt_op | Gt_op | Leq_op | Geq_op | Neq_op) -> popMode;

SPACE_AP:		[ \t\r\n]+ -> skip;

mode AVER_MODE; // Inside of Averaging tag

Averaging_methods:	('R-6' | 'R-3' | 'SUM' | 'CENT' 'E'? 'R'?) -> popMode;

Class_name_AM:		SIMPLE_NAME | '*';

SPACE_AM:		[ \t\r\n]+ -> skip;

mode CLAS_MODE; // Inside of Classification tag

Equ_op_CL:		'=';

Class_name:	(SIMPLE_NAME | '*') -> popMode;

SPACE_CL:		[ \t\r\n]+ -> skip;

mode POTE_MODE; // Inside of Potential tag

Equ_op_PT:		'=';

Potential_types:	('BIHA' 'R'? 'M'? 'O'? 'N'? 'I'? 'C'? | 'LOGN' 'O'? 'R'? 'M'? 'A'? 'L'? | 'SQUA' 'R'? 'E'? '-'? 'W'? 'E'? 'L'? 'L'? | 'SOFT' '-'? 'S'? 'Q'? 'U'? 'A'? 'R'? 'E'? | 'SYMM' 'E'? 'T'? 'R'? 'Y'? | 'HIGH' | '3DPO' | 'HARM' 'O'? 'N'? 'I'? 'C'? | 'MULT' 'I'? 'P'? 'L'? 'E'?) -> popMode;

Class_name_PT:		SIMPLE_NAME | '*';

SPACE_PT:		[ \t\r\n]+ -> skip;

mode ANAL_MODE; // Inside of Noe/Analysis tag

Noe_analysis:		('CURR' 'E'? 'N'? 'T'? | 'TAVE' 'R'? 'A'? 'G'? 'E'? | 'RAVE' 'R'? 'A'? 'G'? 'E'?) -> popMode;

SPACE_NA:		[ \t\r\n]+ -> skip;

mode FLAG_MODE; // Inside of flag statement

Exclude:		'EXCL' 'U'? 'D'? 'E'?;			// Class_name* | Any_class
Include:		'INCL' 'U'? 'D'? 'E'?;			// Class_name*

End_FL:			'END' -> popMode;

Class_name_FL:		SIMPLE_NAME;
Any_class:		'*';

SPACE_FL:		[ \t\r\n]+ -> skip;

mode VECTOR_EXPR_MODE; // vector expression

R_paren_VE:		')' -> popMode;

Equ_op_VE:		'=';
Add_op_VE:		'+';
Sub_op_VE:		'-';
Mul_op_VE:		'*';
Div_op_VE:		'/';
Exp_op_VE:		('^' | '*' '*');
Comma_VE:		',';

Integer_VE:		DECIMAL;
Real_VE:		(DECIMAL | DEC_DOT_DEC) ('E' ('+' | '-')? DECIMAL)?;

Atom_properties_VE:	('B' | 'BCOM' 'P'? | 'CHAR' 'G'? 'E'? | 'CHEM' 'I'? 'C'? 'A'? 'L'? | 'DX' | 'DY' | 'DZ' | 'FBET' 'A'? | 'HARM' 'O'? 'N'? 'I'? 'C'? 'S'? | 'MASS' | 'NAME' | 'Q' | 'QCOM' 'P'? | 'REFX' | 'REFY' | 'REFZ' | 'RESI' 'D'? 'U'? 'E'? | 'RESN' 'A'? 'M'? 'E'? | 'RMSD' | 'SEGI' 'D'? 'E'? 'N'? 'T'? 'I'? 'F'? 'I'? 'E'? 'R'? | 'STORE1' | 'STORE2' | 'STORE3' | 'STORE4' | 'STORE5' | 'STORE6' | 'STORE7' | 'STORE8' | 'STORE9' | 'PSEU' 'D'? 'O'? | 'VX' | 'VY' | 'VZ' | 'X' | 'XCOM' 'P'? | 'Y' | 'YCOM' 'P'? | 'Z' | 'ZCOM' 'P'?);

Abs_VE:			'ABS' -> pushMode(VECTOR_FUNC_MODE);
Acos_VE:		'ACOS' -> pushMode(VECTOR_FUNC_MODE);
Asin_VE:		'ASIN' -> pushMode(VECTOR_FUNC_MODE);
Cos_VE:			'COS' -> pushMode(VECTOR_FUNC_MODE);
Decode_VE:		'DECODE' -> pushMode(VECTOR_FUNC_MODE);
Encode_VE:		'ENCODE' -> pushMode(VECTOR_FUNC_MODE);
Exp_VE:			'EXP' -> pushMode(VECTOR_FUNC_MODE);
Gauss_VE:		'GAUSS' -> pushMode(VECTOR_FUNC_MODE);
Heavy_VE:		'HEAVY' -> pushMode(VECTOR_FUNC_MODE);
Int_VE:			'INT' -> pushMode(VECTOR_FUNC_MODE);
Log10_VE:		'LOG10' -> pushMode(VECTOR_FUNC_MODE);
Log_VE:			'LOG' -> pushMode(VECTOR_FUNC_MODE);
Max_VE:			'MAX' -> pushMode(VECTOR_FUNC_MODE);
Maxw_VE:		'MAXW' -> pushMode(VECTOR_FUNC_MODE);
Min_VE:			'MIN' -> pushMode(VECTOR_FUNC_MODE);
Mod_VE:			'MOD' -> pushMode(VECTOR_FUNC_MODE);
Norm_VE:		'NORM' -> pushMode(VECTOR_FUNC_MODE);
Random_VE:		'RAND' 'O'? 'M'? -> pushMode(VECTOR_FUNC_MODE);
Sign_VE:		'SIGN' -> pushMode(VECTOR_FUNC_MODE);
Sin_VE:			'SIN' -> pushMode(VECTOR_FUNC_MODE);
Sqrt_VE:		'SQRT' -> pushMode(VECTOR_FUNC_MODE);
Tan_VE:			'TAN' -> pushMode(VECTOR_FUNC_MODE);

Symbol_name_VE:		SYMBOL_NAME;
Simple_name_VE:		SIMPLE_NAME;
Double_quote_string_VE:	'"' ~["\r\n]* '"';

SPACE_VE:		[ \t\r\n]+ -> skip;

mode VECTOR_FUNC_MODE; // vector function

L_paren_VF:		'(' -> pushMode(VECTOR_EXPR_MODE);

SPACE_VF:		[ \t\r\n]+ -> skip;

mode VECTOR_SHOW_MODE; // vector show

L_paren_VS:		'(';
R_paren_VS:		')' -> popMode;

Average_VS:		'AVE' 'R'? 'A'? 'G'? 'E'?;
Element_VS:		'ELEM' 'E'? 'N'? 'T'?;
Max_VS:			'MAX';
Min_VS:			'MIN';
Norm_VS:		'NORM';
Rms_VS:			'RMS';
Sum_VS:			'SUM';

Atom_properties_VS:	('B' | 'BCOM' 'P'? | 'CHAR' 'G'? 'E'? | 'DX' | 'DY' | 'DZ' | 'FBET' 'A'? | 'HARM' 'O'? 'N'? 'I'? 'C'? 'S'? | 'MASS' | 'Q' | 'QCOM' 'P'? | 'REFX' | 'REFY' | 'REFZ' | 'RMSD' | 'VX' | 'VY' | 'VZ' | 'X' | 'XCOM' 'P'? | 'Y' | 'YCOM' 'P'? | 'Z' | 'ZCOM' 'P'?);

SPACE_VS:		[ \t\r\n]+ -> skip;

mode CTL_FOR_MODE; // control statement for

L_paren_CF:		'(';
R_paren_CF:		')' -> popMode;
In_CF:			'IN';

Integer_CF:		'-'? DECIMAL;
Real_CF:		('+' | '-')? (DECIMAL | DEC_DOT_DEC) ('E' ('+' | '-')? DECIMAL)?;
Symbol_name_CF:		SYMBOL_NAME;
Simple_name_CF:		SIMPLE_NAME;

SPACE_CF:		[ \t\r\n]+ -> skip;
COMMENT_CF:		'{' (COMMENT_CF | .)*? '}' -> channel(HIDDEN);

mode LOOP_LABEL_MODE; // loop label

Simple_name_LL:		SIMPLE_NAME -> mode(DEFAULT_MODE);

SPACE_LL:		[ \t\r\n]+ -> skip;

