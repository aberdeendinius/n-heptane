species(
    label = 'C#CC([O])C=C=C[O](22646)',
    structure = SMILES('C#CC([O])C=C=C[O]'),
    E0 = (380.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,540,610,2055,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5339,'amu*angstrom^2'), symmetry=1, barrier=(35.2674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53465,'amu*angstrom^2'), symmetry=1, barrier=(35.2846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.189814,0.0753009,-7.96191e-05,3.99907e-08,-7.51561e-12,45926,30.5], Tmin=(100,'K'), Tmax=(1475.35,'K')), NASAPolynomial(coeffs=[21.3901,0.00726947,-7.68433e-07,-1.48615e-11,4.77313e-15,40594.8,-78.4905], Tmin=(1475.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C#CC=O(21959)',
    structure = SMILES('C#CC=O'),
    E0 = (84.2941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,346.083,346.148,1254.41,1569.49,1569.49],'cm^-1')),
        HinderedRotor(inertia=(0.289847,'amu*angstrom^2'), symmetry=1, barrier=(24.6472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3153.34,'J/mol'), sigma=(5.09602,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.54 K, Pc=54.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00307,0.020875,-1.61807e-05,6.30357e-09,-1.00066e-12,10174.9,9.76778], Tmin=(100,'K'), Tmax=(1462.04,'K')), NASAPolynomial(coeffs=[7.33456,0.00902438,-4.02236e-06,7.59526e-10,-5.26541e-14,8908.32,-12.7744], Tmin=(1462.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.2941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC1OC1[C]=C[O](25199)',
    structure = SMILES('C#CC1OC1[C]=C[O]'),
    E0 = (403.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.36598,0.0778268,-8.14651e-05,4.03513e-08,-7.09871e-12,48810.9,33.8742], Tmin=(100,'K'), Tmax=(1764.1,'K')), NASAPolynomial(coeffs=[17.667,0.00601892,3.95478e-06,-1.13625e-09,8.54507e-14,46553.9,-56.1169], Tmin=(1764.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C1OC1C=C=C[O](25249)',
    structure = SMILES('[CH]=C1OC1C=C=C[O]'),
    E0 = (402.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459731,0.0557035,2.09634e-06,-6.60711e-08,3.57213e-11,48564.6,24.9995], Tmin=(100,'K'), Tmax=(925.803,'K')), NASAPolynomial(coeffs=[27.3342,-0.00307645,4.44089e-06,-8.685e-10,5.13291e-14,41131.5,-115.852], Tmin=(925.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(methyleneoxirane) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C1OC=C=CC1[O](25148)',
    structure = SMILES('[CH]=C1OC=C=CC1[O]'),
    E0 = (409.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58857,0.0325357,4.57702e-05,-8.47665e-08,3.42525e-11,49358.5,21.6112], Tmin=(100,'K'), Tmax=(996.192,'K')), NASAPolynomial(coeffs=[17.7304,0.0179179,-7.80134e-06,1.66508e-09,-1.31531e-13,43651.7,-68.7033], Tmin=(996.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'H(8)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CC(=O)C=C=C[O](25250)',
    structure = SMILES('C#CC(=O)C=C=C[O]'),
    E0 = (215.726,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,375,552.5,462.5,1710,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55356,'amu*angstrom^2'), symmetry=1, barrier=(35.7194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.55199,'amu*angstrom^2'), symmetry=1, barrier=(35.6834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931116,0.0665228,-7.56651e-05,4.36963e-08,-9.95025e-12,26057.3,24.1888], Tmin=(100,'K'), Tmax=(1072.6,'K')), NASAPolynomial(coeffs=[13.9679,0.0179047,-7.67301e-06,1.43561e-09,-1.00041e-13,23260.7,-39.6194], Tmin=(1072.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])C=C=C=O(25251)',
    structure = SMILES('C#CC([O])C=C=C=O'),
    E0 = (424.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,540,610,2055,3010,987.5,1337.5,450,1655,2120,512.5,787.5,2175,525,1380,1390,370,380,2900,435,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000761403,'amu*angstrom^2'), symmetry=1, barrier=(8.64501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.374218,'amu*angstrom^2'), symmetry=1, barrier=(8.60402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70758,0.0445641,-3.6698e-05,1.07369e-08,4.85891e-13,51152,10.8847], Tmin=(100,'K'), Tmax=(989.165,'K')), NASAPolynomial(coeffs=[13.0022,0.0108232,-3.62728e-06,6.43656e-10,-4.54292e-14,48333.8,-46.4333], Tmin=(989.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C=C[O](8556)',
    structure = SMILES('[CH]=C=C[O]'),
    E0 = (269.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7401,0.0176211,1.62543e-05,-4.58526e-08,2.291e-11,32513.4,12.64], Tmin=(100,'K'), Tmax=(883.628,'K')), NASAPolynomial(coeffs=[13.5244,-0.00323064,4.17673e-06,-9.22668e-10,6.45049e-14,29515.7,-44.2318], Tmin=(883.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[C]#C(5143)',
    structure = SMILES('[C]#C'),
    E0 = (551.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03853,0.0115449,-2.13263e-05,1.81932e-08,-5.41587e-12,66398,5.96675], Tmin=(100,'K'), Tmax=(1076.58,'K')), NASAPolynomial(coeffs=[4.00845,0.00206817,6.04935e-08,-1.17707e-10,1.2928e-14,66529.5,2.79656], Tmin=(1076.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.936,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = '[O]C=C=CC=O(22476)',
    structure = SMILES('[O]C=C=CC=O'),
    E0 = (-8.09007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.29808,'amu*angstrom^2'), symmetry=1, barrier=(29.8455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63646,0.0447092,-3.72062e-05,1.25679e-08,-1.05205e-12,-881.798,19.0558], Tmin=(100,'K'), Tmax=(1155.06,'K')), NASAPolynomial(coeffs=[14.233,0.010485,-4.96566e-06,1.00356e-09,-7.36229e-14,-4418.66,-46.2447], Tmin=(1155.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.09007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C#CC(O)=C[C]=C[O](25252)',
    structure = SMILES('C#CC(O)=C[C]=C[O]'),
    E0 = (240.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.26348,0.091446,-0.000107587,5.71512e-08,-1.10196e-11,29106,28.5258], Tmin=(100,'K'), Tmax=(1505.1,'K')), NASAPolynomial(coeffs=[26.352,-0.00143153,4.3956e-06,-1.05185e-09,7.68145e-14,23000.3,-108.662], Tmin=(1505.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C#CC([O])C=C=[C]O(25253)',
    structure = SMILES('C#CC([O])C=C=[C]O'),
    E0 = (478.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2175,525,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20443,'amu*angstrom^2'), symmetry=1, barrier=(27.6923,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20541,'amu*angstrom^2'), symmetry=1, barrier=(27.7148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20186,'amu*angstrom^2'), symmetry=1, barrier=(27.633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0252639,0.0798777,-9.61282e-05,5.56794e-08,-1.21965e-11,57733.6,31.6169], Tmin=(100,'K'), Tmax=(1209.36,'K')), NASAPolynomial(coeffs=[20.0033,0.00862569,-1.54255e-06,1.15173e-10,-2.5467e-15,53255.4,-67.303], Tmin=(1209.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(C=CJO)"""),
)

species(
    label = 'C#CC(O)[C]=C=C[O](25254)',
    structure = SMILES('C#CC(O)[C]=C=C[O]'),
    E0 = (387.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,540,610,2055,2175,525,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24354,'amu*angstrom^2'), symmetry=1, barrier=(28.5914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24615,'amu*angstrom^2'), symmetry=1, barrier=(28.6514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24312,'amu*angstrom^2'), symmetry=1, barrier=(28.5818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0348633,0.0790844,-9.41445e-05,5.4188e-08,-1.18587e-11,46810.2,30.3495], Tmin=(100,'K'), Tmax=(1183.04,'K')), NASAPolynomial(coeffs=[19.702,0.00960018,-2.25641e-06,2.72873e-10,-1.42829e-14,42365.9,-66.9554], Tmin=(1183.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])C#C[CH]O(25255)',
    structure = SMILES('C#CC([O])C#C[CH]O'),
    E0 = (458.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,750,770,3400,2100,3615,1277.5,1000,2100,2175,2250,500,525,550,1380,1390,370,380,2900,435,207.604,207.641,1256.69],'cm^-1')),
        HinderedRotor(inertia=(3.79138,'amu*angstrom^2'), symmetry=1, barrier=(115.962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.79136,'amu*angstrom^2'), symmetry=1, barrier=(115.962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.79055,'amu*angstrom^2'), symmetry=1, barrier=(115.962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.78986,'amu*angstrom^2'), symmetry=1, barrier=(115.962,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.347751,0.0791667,-0.000107433,7.39232e-08,-1.93176e-11,55236.3,30.7705], Tmin=(100,'K'), Tmax=(1059.5,'K')), NASAPolynomial(coeffs=[14.8093,0.0152632,-3.78531e-06,4.15334e-10,-1.65538e-14,52694.2,-37.3688], Tmin=(1059.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CCsJOH)"""),
)

species(
    label = '[C]#CC(O)C=C=C[O](25256)',
    structure = SMILES('[C]#CC(O)C=C=C[O]'),
    E0 = (487.259,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21553,'amu*angstrom^2'), symmetry=1, barrier=(27.9473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2162,'amu*angstrom^2'), symmetry=1, barrier=(27.9628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20788,'amu*angstrom^2'), symmetry=1, barrier=(27.7716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0478292,0.0785687,-9.3998e-05,5.47552e-08,-1.20701e-11,58753.2,30.5691], Tmin=(100,'K'), Tmax=(1214.44,'K')), NASAPolynomial(coeffs=[19.1469,0.00977301,-1.75195e-06,1.23644e-10,-1.8493e-15,54548.6,-63.4954], Tmin=(1214.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Acetyl)"""),
)

species(
    label = 'C#CC(O)[CH][C]=C=O(25257)',
    structure = SMILES('C#CC(O)[CH][C]=C=O'),
    E0 = (290.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.229844,0.0786631,-9.25634e-05,5.39845e-08,-1.21869e-11,35124.2,27.1008], Tmin=(100,'K'), Tmax=(1093.4,'K')), NASAPolynomial(coeffs=[17.5559,0.0152787,-5.60809e-06,9.65857e-10,-6.44177e-14,31335.3,-58.0342], Tmin=(1093.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CCJCO) + radical(CCCJ=C=O)"""),
)

species(
    label = 'C#CC([O])=C[C]=CO(25258)',
    structure = SMILES('C#CC([O])=C[C]=CO'),
    E0 = (236.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.17452,0.093013,-0.000114217,6.34665e-08,-1.27967e-11,28659.9,28.2874], Tmin=(100,'K'), Tmax=(1436.93,'K')), NASAPolynomial(coeffs=[25.8423,-0.00162043,4.84912e-06,-1.18303e-09,8.80678e-14,22901.3,-104.868], Tmin=(1436.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=CJC=C)"""),
)

species(
    label = '[C]#CC([O])C=C=CO(25259)',
    structure = SMILES('[C]#CC([O])C=C=CO'),
    E0 = (576.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3615,1277.5,1000,2175,525,1380,1390,370,380,2900,435,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29642,'amu*angstrom^2'), symmetry=1, barrier=(29.8073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28828,'amu*angstrom^2'), symmetry=1, barrier=(29.6201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29398,'amu*angstrom^2'), symmetry=1, barrier=(29.7511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.459053,0.084646,-0.000102629,5.84065e-08,-1.2297e-11,69467.9,30.7826], Tmin=(100,'K'), Tmax=(1334.48,'K')), NASAPolynomial(coeffs=[22.1589,0.00447422,1.39777e-06,-5.11831e-10,4.26057e-14,64533.2,-80.7333], Tmin=(1334.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Acetyl)"""),
)

species(
    label = '[O]C=[C]C=C[O](23191)',
    structure = SMILES('[O]C=[C]C=C[O]'),
    E0 = (158.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62308,'amu*angstrom^2'), symmetry=1, barrier=(37.3178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34928,0.0410469,6.39572e-06,-6.06125e-08,3.33202e-11,19145.6,19.648], Tmin=(100,'K'), Tmax=(893.363,'K')), NASAPolynomial(coeffs=[22.6723,-0.00687368,7.01377e-06,-1.49144e-09,1.02026e-13,13438.2,-91.4392], Tmin=(893.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJC=C) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])=C[C]=C[O](25260)',
    structure = SMILES('C#CC([O])=C[C]=C[O]'),
    E0 = (378.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.60019,'amu*angstrom^2'), symmetry=1, barrier=(36.7916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.60117,'amu*angstrom^2'), symmetry=1, barrier=(36.814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.517494,0.0830948,-0.000100228,5.55071e-08,-1.12818e-11,45646.4,27.3428], Tmin=(100,'K'), Tmax=(1397.4,'K')), NASAPolynomial(coeffs=[23.1998,0.00124401,2.61852e-06,-7.07872e-10,5.429e-14,40380.9,-90.1385], Tmin=(1397.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])[C]=C=C[O](25261)',
    structure = SMILES('C#CC([O])[C]=C=C[O]'),
    E0 = (618.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,540,610,2055,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,1685,370,2175,525,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.52031,'amu*angstrom^2'), symmetry=1, barrier=(34.9549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51751,'amu*angstrom^2'), symmetry=1, barrier=(34.8906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.167356,0.0754383,-8.9409e-05,5.06072e-08,-1.08461e-11,74512.1,29.6826], Tmin=(100,'K'), Tmax=(1220.13,'K')), NASAPolynomial(coeffs=[19.7276,0.00765728,-1.58593e-06,1.65733e-10,-7.68937e-15,70011,-67.4605], Tmin=(1220.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[CH][C]=C=O(25262)',
    structure = SMILES('C#CC([O])[CH][C]=C=O'),
    E0 = (521.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,1380,1390,370,380,2900,435,1685,370,2175,525,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.88489,'amu*angstrom^2'), symmetry=1, barrier=(43.3374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88717,'amu*angstrom^2'), symmetry=1, barrier=(43.3897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8873,'amu*angstrom^2'), symmetry=1, barrier=(43.3927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.369744,0.0749384,-8.75886e-05,5.01334e-08,-1.10728e-11,62825.7,26.4068], Tmin=(100,'K'), Tmax=(1118.42,'K')), NASAPolynomial(coeffs=[17.6062,0.0132935,-4.91296e-06,8.52841e-10,-5.73328e-14,58970.1,-58.6785], Tmin=(1118.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=C=O) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C=C=C[O](25263)',
    structure = SMILES('[C]#CC([O])C=C=C[O]'),
    E0 = (717.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.60455,'amu*angstrom^2'), symmetry=1, barrier=(36.8918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61018,'amu*angstrom^2'), symmetry=1, barrier=(37.0213,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.181819,0.0749068,-8.92162e-05,5.11276e-08,-1.10433e-11,86455,29.8968], Tmin=(100,'K'), Tmax=(1247.49,'K')), NASAPolynomial(coeffs=[19.1255,0.00790338,-1.12118e-06,2.54605e-11,4.02676e-15,82215.8,-63.7317], Tmin=(1247.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(717.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])C=C1[CH]O1(25264)',
    structure = SMILES('C#CC([O])C=C1[CH]O1'),
    E0 = (406.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.833774,0.0471031,2.16972e-05,-8.32248e-08,4.13867e-11,49044.7,25.4011], Tmin=(100,'K'), Tmax=(919.183,'K')), NASAPolynomial(coeffs=[25.3862,-0.000929316,4.10606e-06,-8.57748e-10,5.23409e-14,42046.5,-104.496], Tmin=(919.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(methyleneoxirane) + radical(C=CCJO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=[C]C1C=C(C=O)O1(25147)',
    structure = SMILES('[CH]=[C]C1C=C(C=O)O1'),
    E0 = (426.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.665795,0.0593155,-3.19276e-05,-1.00213e-08,9.65092e-12,51420.5,25.0647], Tmin=(100,'K'), Tmax=(1020.66,'K')), NASAPolynomial(coeffs=[20.5324,0.0112767,-5.15173e-06,1.11409e-09,-8.79486e-14,45811.9,-78.7954], Tmin=(1020.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=CC1[C]=CO1(25125)',
    structure = SMILES('[O]C=C=CC1[C]=CO1'),
    E0 = (373.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.634336,0.0476448,3.11551e-05,-9.88451e-08,4.78507e-11,45108,25.0162], Tmin=(100,'K'), Tmax=(926.512,'K')), NASAPolynomial(coeffs=[28.4972,-0.00456451,5.45739e-06,-1.044e-09,6.08666e-14,37022.8,-123.049], Tmin=(926.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclobutene) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])C1[C]=CO1(25265)',
    structure = SMILES('C#CC([O])C1[C]=CO1'),
    E0 = (504.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615777,0.0501421,2.10461e-05,-8.83768e-08,4.47001e-11,60770.2,26.0536], Tmin=(100,'K'), Tmax=(915.291,'K')), NASAPolynomial(coeffs=[27.7511,-0.00458239,6.071e-06,-1.23951e-09,7.84389e-14,53127.8,-117.07], Tmin=(915.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC1C=[C]C([O])O1(25266)',
    structure = SMILES('C#CC1C=[C]C([O])O1'),
    E0 = (400.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21158,0.0535814,-4.04976e-05,1.52115e-08,-2.28461e-12,48326,24.8679], Tmin=(100,'K'), Tmax=(1572.94,'K')), NASAPolynomial(coeffs=[14.5965,0.0195437,-8.03847e-06,1.45425e-09,-9.80823e-14,44115.3,-45.7692], Tmin=(1572.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + ring(25dihydrofuran) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[O]C1[C]=COC=C=C1(25185)',
    structure = SMILES('[O]C1[C]=COC=C=C1'),
    E0 = (487.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38302,0.0422373,9.80236e-06,-4.81247e-08,2.24607e-11,58688,20.7611], Tmin=(100,'K'), Tmax=(975.204,'K')), NASAPolynomial(coeffs=[16.9168,0.015359,-5.51531e-06,1.0806e-09,-8.30382e-14,53906.6,-62.7712], Tmin=(975.204,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(=O)C=C=CO(25267)',
    structure = SMILES('C#CC(=O)C=C=CO'),
    E0 = (74.2634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.268417,0.0765154,-8.99118e-05,5.19251e-08,-1.1538e-11,9070.99,25.1527], Tmin=(100,'K'), Tmax=(1114.24,'K')), NASAPolynomial(coeffs=[17.9834,0.0129214,-4.30228e-06,7.04407e-10,-4.58156e-14,5123.18,-62.2284], Tmin=(1114.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.2634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC(O)C=C=C=O(25268)',
    structure = SMILES('C#CC(O)C=C=C=O'),
    E0 = (194.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5435,0.0485743,-4.2668e-05,1.58624e-08,-1.16097e-12,23451.6,11.6653], Tmin=(100,'K'), Tmax=(991.513,'K')), NASAPolynomial(coeffs=[13.0741,0.0126038,-4.20569e-06,7.29322e-10,-5.02594e-14,20646.6,-46.4791], Tmin=(991.513,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH)"""),
)

species(
    label = '[CH]=C=C([O])C[C]=C[O](25269)',
    structure = SMILES('[CH]=C=C([O])C[C]=C[O]'),
    E0 = (477.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.967854,'amu*angstrom^2'), symmetry=1, barrier=(22.2529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966398,'amu*angstrom^2'), symmetry=1, barrier=(22.2194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.147834,0.0783578,-9.07721e-05,5.01619e-08,-1.03385e-11,57598,32.4937], Tmin=(100,'K'), Tmax=(1347.66,'K')), NASAPolynomial(coeffs=[20.6651,0.00666782,5.73165e-08,-2.29152e-10,2.21e-14,52888.6,-70.7864], Tmin=(1347.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(C=C=CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CC([O])[C]C=C[O](25270)',
    structure = SMILES('C#CC([O])[C]C=C[O]'),
    E0 = (668.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49182,'amu*angstrom^2'), symmetry=1, barrier=(34.2998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49675,'amu*angstrom^2'), symmetry=1, barrier=(34.4133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49785,'amu*angstrom^2'), symmetry=1, barrier=(34.4385,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.601166,0.079377,-8.54283e-05,4.2776e-08,-7.89192e-12,80537,33.433], Tmin=(100,'K'), Tmax=(1549.05,'K')), NASAPolynomial(coeffs=[22.9911,0.00423989,1.09613e-06,-3.86362e-10,3.0194e-14,74933.5,-85.2053], Tmin=(1549.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CCJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = 'C#CC([O])[CH]C=[C][O](25271)',
    structure = SMILES('C#CC([O])[CH]C=[C][O]'),
    E0 = (576.535,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,1685,370,2175,525,222.211,222.212,222.212,222.213],'cm^-1')),
        HinderedRotor(inertia=(1.16593,'amu*angstrom^2'), symmetry=1, barrier=(40.8537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16592,'amu*angstrom^2'), symmetry=1, barrier=(40.8537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16592,'amu*angstrom^2'), symmetry=1, barrier=(40.8537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.100074,0.0767397,-8.23775e-05,4.24941e-08,-8.29995e-12,69500.2,32.0498], Tmin=(100,'K'), Tmax=(1371.57,'K')), NASAPolynomial(coeffs=[20.937,0.00934367,-2.0603e-06,2.4203e-10,-1.28929e-14,64298,-74.016], Tmin=(1371.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.535,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CJO) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C=C([O])[CH]C=C[O](25272)',
    structure = SMILES('[CH]=C=C([O])[CH]C=C[O]'),
    E0 = (356.638,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63025,'amu*angstrom^2'), symmetry=1, barrier=(37.4826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62972,'amu*angstrom^2'), symmetry=1, barrier=(37.4706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.796488,0.0830031,-9.14812e-05,4.67779e-08,-8.77697e-12,43085.7,31.5589], Tmin=(100,'K'), Tmax=(1533.02,'K')), NASAPolynomial(coeffs=[23.6104,0.00346355,1.85975e-06,-5.60308e-10,4.30321e-14,37465.6,-90.5411], Tmin=(1533.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ) + radical(C=CCJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C=C([O])C=[C]C[O](25273)',
    structure = SMILES('[CH]=C=C([O])C=[C]C[O]'),
    E0 = (587.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.413361,'amu*angstrom^2'), symmetry=1, barrier=(9.50399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413459,'amu*angstrom^2'), symmetry=1, barrier=(9.50624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.48482,0.0845976,-0.000136714,1.15386e-07,-3.72356e-11,70752.5,29.7378], Tmin=(100,'K'), Tmax=(903.431,'K')), NASAPolynomial(coeffs=[9.98807,0.0256749,-1.09113e-05,1.9126e-09,-1.23295e-13,69722.8,-11.3398], Tmin=(903.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=C=CJ) + radical(Cds_S)"""),
)

species(
    label = '[C]#CC([O])C[C]=C[O](25274)',
    structure = SMILES('[C]#CC([O])C[C]=C[O]'),
    E0 = (794.86,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,347.598,347.6,347.6,347.6,347.601],'cm^-1')),
        HinderedRotor(inertia=(0.902833,'amu*angstrom^2'), symmetry=1, barrier=(77.4091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238512,'amu*angstrom^2'), symmetry=1, barrier=(20.4501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.238512,'amu*angstrom^2'), symmetry=1, barrier=(20.4501,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0218937,0.0778529,-9.11867e-05,5.19577e-08,-1.11744e-11,95751.1,32.4913], Tmin=(100,'K'), Tmax=(1258.95,'K')), NASAPolynomial(coeffs=[19.3395,0.00951814,-1.47768e-06,6.29465e-11,2.47783e-15,91438.5,-62.9633], Tmin=(1258.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(794.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(Acetyl) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC([O])[CH]C=C[O](25275)',
    structure = SMILES('[C]#CC([O])[CH]C=C[O]'),
    E0 = (673.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.84555,'amu*angstrom^2'), symmetry=1, barrier=(42.4328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8449,'amu*angstrom^2'), symmetry=1, barrier=(42.4178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.84604,'amu*angstrom^2'), symmetry=1, barrier=(42.4442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.578956,0.0819932,-9.03775e-05,4.69144e-08,-9.02386e-12,81236.5,31.3807], Tmin=(100,'K'), Tmax=(1469.7,'K')), NASAPolynomial(coeffs=[22.6896,0.0057621,5.9366e-07,-3.24024e-10,2.75945e-14,75790.5,-85.0955], Tmin=(1469.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=COJ) + radical(CC(C)OJ) + radical(Acetyl) + radical(C=CCJCO)"""),
)

species(
    label = '[C]#CC([O])C=[C]C[O](25276)',
    structure = SMILES('[C]#CC([O])C=[C]C[O]'),
    E0 = (918.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,180,180,180,1606.7,1608.59],'cm^-1')),
        HinderedRotor(inertia=(0.231132,'amu*angstrom^2'), symmetry=1, barrier=(5.31417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231167,'amu*angstrom^2'), symmetry=1, barrier=(5.31499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.57204,'amu*angstrom^2'), symmetry=1, barrier=(62.1083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.868878,0.0790457,-0.000133265,1.22237e-07,-4.26237e-11,110604,32.2015], Tmin=(100,'K'), Tmax=(881.579,'K')), NASAPolynomial(coeffs=[5.80288,0.0331795,-1.52743e-05,2.79911e-09,-1.86129e-13,110646,14.1944], Tmin=(881.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(918.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(CCOJ) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C=C1[CH]C([O])[C]=C1(25277)',
    structure = SMILES('[O]C=C1[CH]C([O])[C]=C1'),
    E0 = (414.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2950,3150,900,1000,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39286,0.0269324,8.5038e-05,-1.47638e-07,6.30536e-11,50014.5,25.5168], Tmin=(100,'K'), Tmax=(935.42,'K')), NASAPolynomial(coeffs=[25.7925,0.000269776,3.23827e-06,-5.70939e-10,2.39798e-14,42051.4,-108.733], Tmin=(935.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(3-Methylenecyclopentene) + radical(C=CCJC(O)C=C) + radical(cyclopentene-vinyl) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'O(T)(63)',
    structure = SMILES('[O]'),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29230.2,4.09104], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29230.2,4.09104], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CC=C[C]=C[O](22635)',
    structure = SMILES('C#CC=C[C]=C[O]'),
    E0 = (454.145,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62614,'amu*angstrom^2'), symmetry=1, barrier=(37.3881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62728,'amu*angstrom^2'), symmetry=1, barrier=(37.4144,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831519,0.0562972,-2.5861e-05,-2.38426e-08,1.8337e-11,54747.6,21.7645], Tmin=(100,'K'), Tmax=(911.449,'K')), NASAPolynomial(coeffs=[20.957,0.00392216,1.17357e-06,-3.44429e-10,2.26331e-14,49585.7,-81.6544], Tmin=(911.449,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C=CC([O])C#C(22633)',
    structure = SMILES('[CH]=C=CC([O])C#C'),
    E0 = (602.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,540,610,2055,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2175,525,180],'cm^-1')),
        HinderedRotor(inertia=(1.68171,'amu*angstrom^2'), symmetry=1, barrier=(38.6658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68241,'amu*angstrom^2'), symmetry=1, barrier=(38.6818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09961,0.0599246,-6.02971e-05,2.72634e-08,-3.29184e-12,72546,25.2234], Tmin=(100,'K'), Tmax=(898.987,'K')), NASAPolynomial(coeffs=[14.0593,0.0143154,-4.30948e-06,6.5984e-10,-4.14889e-14,69728.7,-38.6287], Tmin=(898.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Cdd-CdsCds) + group(Ct-CtH) + radical(CC(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C1C(C=O)=CC1[O](25066)',
    structure = SMILES('[CH]=C1C(C=O)=CC1[O]'),
    E0 = (405.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752662,0.0551801,-1.58823e-05,-2.96103e-08,1.74329e-11,48889.1,24.0029], Tmin=(100,'K'), Tmax=(986.81,'K')), NASAPolynomial(coeffs=[21.3006,0.00952956,-3.7057e-06,8.16256e-10,-6.78155e-14,43001,-84.1422], Tmin=(986.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])C#CC=O(25278)',
    structure = SMILES('C#CC([O])C#CC=O'),
    E0 = (343.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2100,2175,2250,500,525,550,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.80068,'amu*angstrom^2'), symmetry=1, barrier=(41.4012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80928,'amu*angstrom^2'), symmetry=1, barrier=(41.599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.80417,'amu*angstrom^2'), symmetry=1, barrier=(41.4815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.923409,0.0650628,-7.43186e-05,4.37208e-08,-1.00814e-11,41392.5,26.8123], Tmin=(100,'K'), Tmax=(1065.64,'K')), NASAPolynomial(coeffs=[13.8346,0.0165982,-6.09849e-06,1.04151e-09,-6.8666e-14,38640.8,-36.2976], Tmin=(1065.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[C]=CC=O(25279)',
    structure = SMILES('C#CC([O])[C]=CC=O'),
    E0 = (421.581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,2175,525,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.917075,'amu*angstrom^2'), symmetry=1, barrier=(21.0854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.916951,'amu*angstrom^2'), symmetry=1, barrier=(21.0825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91685,'amu*angstrom^2'), symmetry=1, barrier=(21.0802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0283,0.0708232,-8.87366e-05,6.10218e-08,-1.7261e-11,50806.7,27.5608], Tmin=(100,'K'), Tmax=(852.555,'K')), NASAPolynomial(coeffs=[10.0388,0.0285471,-1.43538e-05,2.85609e-09,-2.04411e-13,49270.4,-14.4722], Tmin=(852.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[CH]C=C=O(25280)',
    structure = SMILES('C#CC([O])[CH]C=C=O'),
    E0 = (318.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.363621,0.074129,-8.1071e-05,4.40767e-08,-9.31074e-12,38497.9,26.7201], Tmin=(100,'K'), Tmax=(1164.56,'K')), NASAPolynomial(coeffs=[17.3618,0.0157445,-5.8703e-06,1.02767e-09,-6.93565e-14,34538.7,-57.8763], Tmin=(1164.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CCJCO) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])=CC=C[O](25281)',
    structure = SMILES('C#CC([O])=CC=C[O]'),
    E0 = (179.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.948256,0.0822077,-8.94285e-05,4.44229e-08,-8.03931e-12,21737.2,29.3159], Tmin=(100,'K'), Tmax=(1608.07,'K')), NASAPolynomial(coeffs=[24.3124,0.00144969,2.62091e-06,-6.69693e-10,4.86146e-14,15930.4,-97.3462], Tmin=(1608.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=COJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C=CC=O(25282)',
    structure = SMILES('[C]#CC([O])C=CC=O'),
    E0 = (520.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,2175,525,2995,3025,975,1000,1300,1375,400,500,1630,1680,291.518,291.527,291.53],'cm^-1')),
        HinderedRotor(inertia=(0.225916,'amu*angstrom^2'), symmetry=1, barrier=(13.6258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225895,'amu*angstrom^2'), symmetry=1, barrier=(13.6259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716224,'amu*angstrom^2'), symmetry=1, barrier=(43.1929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.991446,0.0709246,-9.08705e-05,6.46863e-08,-1.88454e-11,62751.9,27.957], Tmin=(100,'K'), Tmax=(831.6,'K')), NASAPolynomial(coeffs=[9.91741,0.0279919,-1.34329e-05,2.60898e-09,-1.83887e-13,61267.2,-13.4599], Tmin=(831.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(520.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C1[C]=CC(C=O)=C1(25094)',
    structure = SMILES('[O]C1[C]=CC(C=O)=C1'),
    E0 = (331.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13647,0.0490379,-1.10067e-05,-2.48062e-08,1.32015e-11,40009.5,23.3329], Tmin=(100,'K'), Tmax=(1033.5,'K')), NASAPolynomial(coeffs=[17.6442,0.0158228,-7.32088e-06,1.53551e-09,-1.17564e-13,34959.1,-64.7772], Tmin=(1033.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopentadiene) + radical(1,3-cyclopentadiene-vinyl-1) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=[C]C1C=C=COO1(25214)',
    structure = SMILES('[CH]=[C]C1C=C=COO1'),
    E0 = (701.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.909029,0.11813,-0.000180299,1.36136e-07,-3.84464e-11,84546.9,15.7001], Tmin=(100,'K'), Tmax=(662.628,'K')), NASAPolynomial(coeffs=[15.622,0.0321699,-1.7019e-05,3.36069e-09,-2.36373e-13,82052.5,-59.5405], Tmin=(662.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = 'C#CC(=O)C=CC=O(25283)',
    structure = SMILES('C#CC(=O)C=CC=O'),
    E0 = (18.9893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15887,0.0698004,-0.000104872,9.68217e-08,-3.67439e-11,2379.11,24.3123], Tmin=(100,'K'), Tmax=(759.232,'K')), NASAPolynomial(coeffs=[5.49876,0.0366904,-1.92155e-05,3.83447e-09,-2.72461e-13,2015.41,6.5149], Tmin=(759.232,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.9893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'C#CC1C=C(C=O)O1(25170)',
    structure = SMILES('C#CC1C=C(C=O)O1'),
    E0 = (107.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.764186,0.053925,-9.20228e-06,-3.9698e-08,2.21183e-11,13054.9,21.4488], Tmin=(100,'K'), Tmax=(961.894,'K')), NASAPolynomial(coeffs=[21.9519,0.00756437,-2.00843e-06,4.36876e-10,-3.99684e-14,7047.48,-89.9855], Tmin=(961.894,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=O(373)',
    structure = SMILES('[CH]=O'),
    E0 = (33.3613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([880.033,2045.55,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.35602,-0.00347087,1.25664e-05,-9.99488e-09,2.27888e-12,3995.77,2.75113], Tmin=(100,'K'), Tmax=(1565.72,'K')), NASAPolynomial(coeffs=[4.61867,0.00504456,-4.39241e-06,9.73282e-10,-7.07436e-14,2787.5,-2.22961], Tmin=(1565.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(HCdsJO)"""),
)

species(
    label = '[C]=CC([O])C#C(23423)',
    structure = SMILES('[C]=CC([O])C#C'),
    E0 = (865.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,247.16,248.276],'cm^-1')),
        HinderedRotor(inertia=(0.904632,'amu*angstrom^2'), symmetry=1, barrier=(40.45,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.497215,'amu*angstrom^2'), symmetry=1, barrier=(22.6292,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73191,0.0491992,-5.68503e-05,3.40177e-08,-8.03344e-12,104157,21.8463], Tmin=(100,'K'), Tmax=(1036.5,'K')), NASAPolynomial(coeffs=[10.9264,0.0137164,-5.50057e-06,9.90213e-10,-6.73574e-14,102251,-22.8416], Tmin=(1036.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(865.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (380.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (464.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (464.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (557.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (449.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (639.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (380.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (602.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (505.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (656.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (529.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (773.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (571.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (471.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (471.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (658.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (539.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (710.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (589.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (830.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (733.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (929.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (567.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (506.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (622.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (506.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (450.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (511.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (405.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (405.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (505.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (690.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (599.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (445.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (626.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (812.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (698.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (943.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (565.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (861.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (1009.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (434.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (564.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (618.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (543.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (557.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (594.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (436.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (709.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (469.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (388.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (933.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC=O(21959)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC1OC1[C]=C[O](25199)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(83.6851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[CH]=C1OC1C=C=C[O](25249)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_T;triplebond_intra_H;radadd_intra] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[CH]=C1OC=C=CC1[O](25148)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.05e+09,'s^-1'), n=0.155, Ea=(177.192,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;triplebond_intra_H;radadd_intra] for rate rule [R7_SMMS;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C#CC(=O)C=C=C[O](25250)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;HJ] for rate rule [CO-DeDe_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', 'C#CC([O])C=C=C=O(25251)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=C[O](8556)', 'C#CC=O(21959)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.99e+06,'cm^3/(mol*s)'), n=2.12, Ea=(26.3073,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CtH_O;YJ] for rate rule [CO-CtH_O;CdsJ=Cdd]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 21.2 to 26.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[C]#C(5143)', '[O]C=C=CC=O(22476)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CdH_O;CtJ_Ct]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC(O)=C[C]=C[O](25252)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C#CC([O])C=C=[C]O(25253)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CC(O)[C]=C=C[O](25254)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([O])C#C[CH]O(25255)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[C]#CC(O)C=C=C[O](25256)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;O_H_out] for rate rule [R4H_TSS;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC(O)[CH][C]=C=O(25257)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC([O])=C[C]=CO(25258)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[C]#CC([O])C=C=CO(25259)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.21847e+06,'s^-1'), n=1.22418, Ea=(82.4275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7H;Y_rad_out;XH_out] for rate rule [R7H;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C=C[O](8556)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(81155.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]#C(5143)', '[O]C=[C]C=C[O](23191)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.34536e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_rad/Ct;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', 'C#CC([O])=C[C]=C[O](25260)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', 'C#CC([O])[C]=C=C[O](25261)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', 'C#CC([O])[CH][C]=C=O(25262)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', '[C]#CC([O])C=C=C[O](25263)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC([O])C=C1[CH]O1(25264)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[CH]=[C]C1C=C(C=O)O1(25147)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[O]C=C=CC1[C]=CO1(25125)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC([O])C1[C]=CO1(25265)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC1C=[C]C([O])O1(25266)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(69.6427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[O]C1[C]=COC=C=C1(25185)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC(=O)C=C=CO(25267)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC(O)C=C=C=O(25268)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C=C([O])C[C]=C[O](25269)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([O])[C]C=C[O](25270)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([O])[CH]C=[C][O](25271)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C=C([O])[CH]C=C[O](25272)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=C=C([O])C=[C]C[O](25273)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#CC([O])C[C]=C[O](25274)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[C]#CC([O])[CH]C=C[O](25275)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[C]#CC([O])C=[C]C[O](25276)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]C=C1[CH]C([O])[C]=C1(25277)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O(T)(63)', 'C#CC=C[C]=C[O](22635)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/TwoDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['O(T)(63)', '[CH]=C=CC([O])C#C(22633)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[CH]=C1C(C=O)=CC1[O](25066)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(8)', 'C#CC([O])C#CC=O(25278)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C#CC([O])[C]=CC=O(25279)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC([O])[CH]C=C=O(25280)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC([O])=CC=C[O](25281)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(8.2826e+06,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[C]#CC([O])C=CC=O(25282)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(272058,'s^-1'), n=1.7475, Ea=(73.8225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cd_H_out_single] for rate rule [R5H_TSSD;Ct_rad_out;Cd_H_out_singleDe]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[O]C1[C]=CC(C=O)=C1(25094)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(8.72e+09,'s^-1'), n=0.186, Ea=(55.7727,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 831 used for R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleDe
Exact match found for rate rule [R5_DS_T;triplebond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['[CH]=[C]C1C=C=COO1(25214)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(329.213,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC(=O)C=CC=O(25283)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C#CC([O])C=C=C[O](22646)'],
    products = ['C#CC1C=C(C=O)O1(25170)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=O(373)', '[C]=CC([O])C#C(23423)'],
    products = ['C#CC([O])C=C=C[O](22646)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4949',
    isomers = [
        'C#CC([O])C=C=C[O](22646)',
    ],
    reactants = [
        ('C#CC=O(21959)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4949',
    Tmin = (1200,'K'),
    Tmax = (1500,'K'),
    Tcount = 10,
    Tlist = ([1201.48,1213.22,1236.21,1269.31,1310.55,1356.92,1404.16,1447.02,1479.84,1497.7],'K'),
    Pmin = (1,'atm'),
    Pmax = (10,'atm'),
    Pcount = 10,
    Plist = ([1.02771,1.14872,1.41959,1.89986,2.67608,3.83649,5.40396,7.23219,8.93758,9.98989],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

