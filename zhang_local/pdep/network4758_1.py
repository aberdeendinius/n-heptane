species(
    label = '[CH]=C=COO[C]=C(22461)',
    structure = SMILES('[CH]=C=COO[C]=C'),
    E0 = (597.945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.03086,'amu*angstrom^2'), symmetry=1, barrier=(23.7015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03051,'amu*angstrom^2'), symmetry=1, barrier=(23.6934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02858,'amu*angstrom^2'), symmetry=1, barrier=(23.649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15813,0.0619509,-7.14597e-05,4.30466e-08,-1.02667e-11,72019,28.8467], Tmin=(100,'K'), Tmax=(1025.19,'K')), NASAPolynomial(coeffs=[12.4005,0.0180862,-7.27925e-06,1.31093e-09,-8.91114e-14,69713.9,-25.6711], Tmin=(1025.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C=O(598)',
    structure = SMILES('C=C=O'),
    E0 = (-59.3981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52746,0.00708371,9.17709e-06,-1.64254e-08,6.71115e-12,-7123.94,5.7438], Tmin=(100,'K'), Tmax=(956.683,'K')), NASAPolynomial(coeffs=[5.76495,0.00596559,-1.98486e-06,3.52744e-10,-2.51619e-14,-7929,-6.92178], Tmin=(956.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-59.3981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-(Cdd-O2d)HH)"""),
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
    label = '[CH]=[C]C1OOC1=C(23818)',
    structure = SMILES('[CH]=[C]C1OOC1=C'),
    E0 = (633.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78264,0.0369184,2.96124e-06,-3.02077e-08,1.34827e-11,76310.5,24.6092], Tmin=(100,'K'), Tmax=(1040.92,'K')), NASAPolynomial(coeffs=[13.7799,0.0176863,-8.04641e-06,1.64207e-09,-1.22969e-13,72357.1,-40.7444], Tmin=(1040.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.732,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=C=COOC#C(23819)',
    structure = SMILES('[CH]=C=COOC#C'),
    E0 = (589.717,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,750,770,3400,2100,540,610,2055,350,500,795,815,2175,525],'cm^-1')),
        HinderedRotor(inertia=(1.61677,'amu*angstrom^2'), symmetry=1, barrier=(37.1727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6176,'amu*angstrom^2'), symmetry=1, barrier=(37.1919,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851836,0.0675712,-8.41343e-05,5.17005e-08,-1.2303e-11,71041.5,24.4247], Tmin=(100,'K'), Tmax=(1037.22,'K')), NASAPolynomial(coeffs=[15.0192,0.0129348,-5.12026e-06,9.14469e-10,-6.20644e-14,68102.6,-44.4426], Tmin=(1037.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(589.717,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=C=CJ)"""),
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
    label = '[CH]=C=COOC=[CH](23820)',
    structure = SMILES('[CH]=C=COOC=[CH]'),
    E0 = (605.297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,540,610,2055,350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(1.27599,'amu*angstrom^2'), symmetry=1, barrier=(29.3375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2661,'amu*angstrom^2'), symmetry=1, barrier=(29.1101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.506102,0.0681032,-7.67651e-05,4.20328e-08,-8.7654e-12,72933.7,28.4577], Tmin=(100,'K'), Tmax=(1250.38,'K')), NASAPolynomial(coeffs=[17.9426,0.00916525,-2.27258e-06,2.95591e-10,-1.66214e-14,68820.1,-58.5723], Tmin=(1250.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C#COO[C]=C(23821)',
    structure = SMILES('[CH2]C#COO[C]=C'),
    E0 = (648.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,350,500,795,815,2100,2250,500,550,3000,3100,440,815,1455,1000,327.458],'cm^-1')),
        HinderedRotor(inertia=(0.342168,'amu*angstrom^2'), symmetry=1, barrier=(26.0579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721233,'amu*angstrom^2'), symmetry=1, barrier=(5.48987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.560404,'amu*angstrom^2'), symmetry=1, barrier=(42.667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.560365,'amu*angstrom^2'), symmetry=1, barrier=(42.6631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64127,0.0564459,-6.20127e-05,3.86525e-08,-1.02066e-11,78093.5,27.7263], Tmin=(100,'K'), Tmax=(897.707,'K')), NASAPolynomial(coeffs=[8.14466,0.0274682,-1.35934e-05,2.69495e-09,-1.92879e-13,76925.9,-2.94711], Tmin=(897.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cs-CtHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(Propargyl) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C=[C]OOC=C(23822)',
    structure = SMILES('[CH]=C=[C]OOC=C'),
    E0 = (597.945,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.03086,'amu*angstrom^2'), symmetry=1, barrier=(23.7015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03051,'amu*angstrom^2'), symmetry=1, barrier=(23.6934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02858,'amu*angstrom^2'), symmetry=1, barrier=(23.649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15813,0.0619509,-7.14597e-05,4.30466e-08,-1.02667e-11,72019,28.8467], Tmin=(100,'K'), Tmax=(1025.19,'K')), NASAPolynomial(coeffs=[12.4005,0.0180862,-7.27925e-06,1.31093e-09,-8.91114e-14,69713.9,-25.6711], Tmin=(1025.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.945,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OOC=C=C(23823)',
    structure = SMILES('[CH]=[C]OOC=C=C'),
    E0 = (690.564,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.949476,'amu*angstrom^2'), symmetry=1, barrier=(21.8303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.950021,'amu*angstrom^2'), symmetry=1, barrier=(21.8429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.949466,'amu*angstrom^2'), symmetry=1, barrier=(21.8301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25886,0.0621883,-7.1849e-05,4.33441e-08,-1.04991e-11,83152.9,28.2675], Tmin=(100,'K'), Tmax=(1000.36,'K')), NASAPolynomial(coeffs=[11.6428,0.0206689,-9.5943e-06,1.85718e-09,-1.31477e-13,81075.3,-21.833], Tmin=(1000.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.564,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=O(601)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,672.051,672.102],'cm^-1')),
        HinderedRotor(inertia=(0.000373196,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57974,0.00389613,2.17609e-05,-3.06386e-08,1.18311e-11,19367.5,10.1348], Tmin=(100,'K'), Tmax=(961.532,'K')), NASAPolynomial(coeffs=[6.4326,0.00553733,-1.87382e-06,3.59985e-10,-2.76653e-14,18194.3,-6.76404], Tmin=(961.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH) + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'C=[C]O[O](591)',
    structure = SMILES('C=[C]O[O]'),
    E0 = (349.6,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.227153,'amu*angstrom^2'), symmetry=1, barrier=(5.22269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.10707,0.0218976,-3.08154e-05,2.70444e-08,-9.38388e-12,42077.2,16.3544], Tmin=(100,'K'), Tmax=(876.588,'K')), NASAPolynomial(coeffs=[4.1919,0.0119052,-5.08868e-06,9.16781e-10,-6.09557e-14,42080.7,12.3686], Tmin=(876.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.6,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C=[CH](18734)',
    structure = SMILES('[CH]=C=[CH]'),
    E0 = (491.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,239.877,511.233,1743.98,1746.51,1747.6,1753.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.766,0.0170203,-1.57568e-05,7.95984e-09,-1.4265e-12,59188.9,11.2142], Tmin=(100,'K'), Tmax=(1806.04,'K')), NASAPolynomial(coeffs=[4.81405,0.00509933,2.77647e-07,-2.23082e-10,1.96202e-14,59653.5,3.45727], Tmin=(1806.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C=[C]OO[C]=C(23824)',
    structure = SMILES('[CH]=C=[C]OO[C]=C'),
    E0 = (837.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,350,500,795,815,1670,1700,300,440],'cm^-1')),
        HinderedRotor(inertia=(0.502597,'amu*angstrom^2'), symmetry=1, barrier=(11.5557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.502609,'amu*angstrom^2'), symmetry=1, barrier=(11.556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.70253,'amu*angstrom^2'), symmetry=1, barrier=(62.1364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31974,0.0643124,-0.00010066,8.52395e-08,-2.81743e-11,100842,31.3932], Tmin=(100,'K'), Tmax=(863.071,'K')), NASAPolynomial(coeffs=[8.22725,0.0222697,-1.01603e-05,1.87089e-09,-1.25532e-13,100023,1.24957], Tmin=(863.071,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(837.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OOC=C=[CH](23825)',
    structure = SMILES('[CH]=[C]OOC=C=[CH]'),
    E0 = (845.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3115,3125,620,680,785,800,1600,1700,1685,370,540,610,2055,350,500,795,815,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.03754,'amu*angstrom^2'), symmetry=1, barrier=(23.8552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04443,'amu*angstrom^2'), symmetry=1, barrier=(24.0134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07097,0.0656171,-8.85604e-05,6.09449e-08,-1.63878e-11,101739,29.5637], Tmin=(100,'K'), Tmax=(916.799,'K')), NASAPolynomial(coeffs=[12.5944,0.0153393,-6.29798e-06,1.1252e-09,-7.53623e-14,99626.5,-25.029], Tmin=(916.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OOC1[C]=C1(23826)',
    structure = SMILES('C=[C]OOC1[C]=C1'),
    E0 = (741.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68789,0.0533232,-4.97404e-05,2.40529e-08,-4.8194e-12,89272.6,26.5734], Tmin=(100,'K'), Tmax=(1167.58,'K')), NASAPolynomial(coeffs=[10.3347,0.0237001,-1.16834e-05,2.32302e-09,-1.66639e-13,87253.4,-16.4821], Tmin=(1167.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=CJO) + radical(cyclopropenyl-vinyl)"""),
)

species(
    label = '[CH]=C1[CH]OOC1=C(23827)',
    structure = SMILES('[CH]=C1[CH]OOC1=C'),
    E0 = (424.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.002,0.0214793,6.74938e-05,-1.11878e-07,4.65042e-11,51131.6,20.4622], Tmin=(100,'K'), Tmax=(944.873,'K')), NASAPolynomial(coeffs=[18.6976,0.00772308,-1.03316e-06,2.30046e-10,-2.76569e-14,45435.5,-72.5837], Tmin=(944.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.363,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Cds_P) + radical(C=CCJO)"""),
)

species(
    label = 'C1=COC=1(22275)',
    structure = SMILES('C1=COC=1'),
    E0 = (388.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17882,0.00883829,2.87655e-05,-4.70958e-08,1.95464e-11,46730.5,10.2865], Tmin=(100,'K'), Tmax=(941.322,'K')), NASAPolynomial(coeffs=[9.87744,0.00380291,-5.45265e-07,1.04211e-10,-1.15627e-14,44431.4,-27.139], Tmin=(941.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[CH]=C(C=O)O[C]=C(22463)',
    structure = SMILES('[CH]=C(C=O)O[C]=C'),
    E0 = (335.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,350,440,435,1725,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,298.999,299.005],'cm^-1')),
        HinderedRotor(inertia=(0.220416,'amu*angstrom^2'), symmetry=1, barrier=(13.9836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220417,'amu*angstrom^2'), symmetry=1, barrier=(13.9835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220411,'amu*angstrom^2'), symmetry=1, barrier=(13.9836,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3604.33,'J/mol'), sigma=(5.78855,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.99 K, Pc=42.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42786,0.0620956,-7.86306e-05,5.53213e-08,-1.61605e-11,40427.3,25.7072], Tmin=(100,'K'), Tmax=(823.082,'K')), NASAPolynomial(coeffs=[8.7523,0.0264996,-1.37582e-05,2.77575e-09,-2.00143e-13,39221.6,-8.20273], Tmin=(823.082,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]=C(584)',
    structure = SMILES('[C]=C'),
    E0 = (600.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.94093,-0.00117598,1.80376e-05,-2.01208e-08,6.96659e-12,72197.9,5.25681], Tmin=(100,'K'), Tmax=(976.125,'K')), NASAPolynomial(coeffs=[3.93016,0.00536132,-1.98619e-06,3.69549e-10,-2.66221e-14,71890.7,3.724], Tmin=(976.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C=CO[O](20803)',
    structure = SMILES('[CH]=C=CO[O]'),
    E0 = (404.912,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,540,610,2055,180,180,180,2867.84],'cm^-1')),
        HinderedRotor(inertia=(0.0105329,'amu*angstrom^2'), symmetry=1, barrier=(9.09612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.1,'J/mol'), sigma=(5.7752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.83 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13751,0.0383774,-4.88548e-05,3.09451e-08,-7.37348e-12,48769,18.1869], Tmin=(100,'K'), Tmax=(1167.31,'K')), NASAPolynomial(coeffs=[10.2545,0.00577971,-8.1986e-07,1.20449e-12,5.54427e-15,47199.8,-20.8326], Tmin=(1167.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(ROOJ)"""),
)

species(
    label = '[C]#C[CH]OO[C]=C(23828)',
    structure = SMILES('[C]#C[CH]OO[C]=C'),
    E0 = (959.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,350,500,795,815,2175,525,2950,3100,1380,975,1025,1650,180,1838.32],'cm^-1')),
        HinderedRotor(inertia=(0.16148,'amu*angstrom^2'), symmetry=1, barrier=(3.71274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72545,'amu*angstrom^2'), symmetry=1, barrier=(62.6634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72465,'amu*angstrom^2'), symmetry=1, barrier=(62.6452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7226,'amu*angstrom^2'), symmetry=1, barrier=(62.598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18278,0.0703909,-0.00012004,1.07484e-07,-3.63616e-11,115442,28.3445], Tmin=(100,'K'), Tmax=(899.142,'K')), NASAPolynomial(coeffs=[6.80728,0.0252882,-1.1297e-05,2.01811e-09,-1.31334e-13,115243,6.32142], Tmin=(899.142,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(959.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOOC) + radical(Acetyl) + radical(C=CJO)"""),
)

species(
    label = '[C]#CCOO[C]=C(23829)',
    structure = SMILES('[C]#CCOO[C]=C'),
    E0 = (772.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,350,500,795,815,2175,525,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.143873,'amu*angstrom^2'), symmetry=1, barrier=(3.30792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144667,'amu*angstrom^2'), symmetry=1, barrier=(3.32618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.62675,'amu*angstrom^2'), symmetry=1, barrier=(60.3941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.6268,'amu*angstrom^2'), symmetry=1, barrier=(60.3953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16285,0.0686904,-0.000106416,9.0942e-08,-3.0107e-11,93021.5,28.1211], Tmin=(100,'K'), Tmax=(888.777,'K')), NASAPolynomial(coeffs=[7.46662,0.0263698,-1.1447e-05,2.04649e-09,-1.34468e-13,92451.9,1.55194], Tmin=(888.777,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(772.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO) + radical(Acetyl)"""),
)

species(
    label = '[CH]=[C]OOCC#C(23830)',
    structure = SMILES('[CH]=[C]OOCC#C'),
    E0 = (682.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,750,770,3400,2100,350,500,795,815,2750,2850,1437.5,1250,1305,750,350,1685,370,2175,525],'cm^-1')),
        HinderedRotor(inertia=(0.483318,'amu*angstrom^2'), symmetry=1, barrier=(11.1124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.90066,'amu*angstrom^2'), symmetry=1, barrier=(89.6838,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08946,0.0683856,-9.94763e-05,8.01248e-08,-2.57021e-11,82195.7,28.0463], Tmin=(100,'K'), Tmax=(841.414,'K')), NASAPolynomial(coeffs=[9.18077,0.0243066,-1.08882e-05,2.00554e-09,-1.35418e-13,81032.8,-8.41161], Tmin=(841.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(682.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[C]#C[CH]OOC=C(23831)',
    structure = SMILES('[C]#C[CH]OOC=C'),
    E0 = (719.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,350,500,795,815,2175,525,3010,987.5,1337.5,450,1655,180,1049.53],'cm^-1')),
        HinderedRotor(inertia=(0.758594,'amu*angstrom^2'), symmetry=1, barrier=(17.4416,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84219,'amu*angstrom^2'), symmetry=1, barrier=(65.3476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84264,'amu*angstrom^2'), symmetry=1, barrier=(65.3578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84277,'amu*angstrom^2'), symmetry=1, barrier=(65.3608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15497,0.0661339,-8.24343e-05,5.11119e-08,-1.04835e-11,86613.7,25.3367], Tmin=(100,'K'), Tmax=(702.477,'K')), NASAPolynomial(coeffs=[10.7609,0.0215116,-8.66521e-06,1.51965e-09,-1.00184e-13,85015.5,-19.3832], Tmin=(702.477,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(719.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCsJOOC) + radical(Acetyl)"""),
)

species(
    label = 'C=C1C=[C][CH]OO1(23832)',
    structure = SMILES('C=C1C=[C][CH]OO1'),
    E0 = (390.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18997,0.0196406,6.67451e-05,-9.85643e-08,3.73661e-11,47020.5,16.8361], Tmin=(100,'K'), Tmax=(1003.91,'K')), NASAPolynomial(coeffs=[14.8291,0.0198306,-9.06748e-06,1.93683e-09,-1.51489e-13,41935.5,-56.8766], Tmin=(1003.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = 'C#CC1OOC1=C(23833)',
    structure = SMILES('C#CC1OOC1=C'),
    E0 = (314.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87397,0.0316314,2.52204e-05,-5.91133e-08,2.55409e-11,37945.2,21.0174], Tmin=(100,'K'), Tmax=(970.536,'K')), NASAPolynomial(coeffs=[15.1241,0.0140929,-4.96813e-06,9.79599e-10,-7.61741e-14,33627.3,-51.5055], Tmin=(970.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Cyclobutane)"""),
)

species(
    label = '[CH]OO[C]=C(4541)',
    structure = SMILES('[CH]OO[C]=C'),
    E0 = (658.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,500,795,815,1685,370,180,793.866,793.887],'cm^-1')),
        HinderedRotor(inertia=(0.144155,'amu*angstrom^2'), symmetry=1, barrier=(3.31441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144113,'amu*angstrom^2'), symmetry=1, barrier=(3.31344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144138,'amu*angstrom^2'), symmetry=1, barrier=(3.31402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20537,0.0431336,-6.10799e-05,4.96098e-08,-1.66725e-11,79261.4,20.9229], Tmin=(100,'K'), Tmax=(721.447,'K')), NASAPolynomial(coeffs=[6.68848,0.018281,-9.41505e-06,1.87487e-09,-1.33513e-13,78614.5,0.75756], Tmin=(721.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(658.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CH2_triplet)"""),
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
    E0 = (597.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (719.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (816.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (597.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (710.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (814.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (743.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (881.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (597.945,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (856.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1049.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1056.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (747.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (641.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (680.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (911.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1010.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1170.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (925.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (726.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (888.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (628.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (605.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1244.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['C=C=O(598)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['[CH]=[C]C1OOC1=C(23818)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH]=C=COOC#C(23819)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C=O(598)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.3e+11,'cm^3/(mol*s)'), n=0, Ea=(387.468,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_R;O_rad/OneDe] for rate rule [Od_Cdd;O_rad/OneDe]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 383.3 to 387.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C=COOC=[CH](23820)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C#COO[C]=C(23821)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.65028e+09,'s^-1'), n=1.32317, Ea=(166.12,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] + [R3H;Cd_rad_out_singleNd;XH_out] for rate rule [R3H;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C=[C]OOC=C(23822)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cd_H_out_doubleC] for rate rule [R4H_SSS;Cd_rad_out_double;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]OOC=C=C(23823)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R7HJ_1;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=O(601)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.63841e+06,'m^3/(mol*s)'), n=0.151, Ea=(167.206,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [O_rad/OneDe;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 167.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]O[O](591)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.15767e+07,'m^3/(mol*s)'), n=0.0716491, Ea=(15.4197,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;O_rad/NonDe] for rate rule [Cd_allenic;O_rad/NonDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(8)', '[CH]=C=[C]OO[C]=C(23824)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(8)', '[CH]=[C]OOC=C=[CH](23825)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['C=[C]OOC1[C]=C1(23826)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.2354e+12,'s^-1'), n=-0.1205, Ea=(149.998,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cdsingleH] + [R3;doublebond_intra_CdCdd;radadd_intra_cdsingle] for rate rule [R3;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['[CH]=C1[CH]OOC1=C(23827)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.89256e+09,'s^-1'), n=0.56, Ea=(43.3044,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['[CH2][C]=O(601)', 'C1=COC=1(22275)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO_intra] for rate rule [R3OO;Cd_pri_rad_in;OO_intra]
Euclidian distance = 2.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['[CH]=C(C=O)O[C]=C(22463)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[C]=C(584)', '[CH]=C=CO[O](20803)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(43.5839,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[C]#C[CH]OO[C]=C(23828)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[C]#CCOO[C]=C(23829)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.9263e+09,'s^-1'), n=1.08337, Ea=(153.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;Cs_H_out] for rate rule [R3H_TS;Ct_rad_out;Cs_H_out_OOH/H]
Euclidian distance = 2.44948974278
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]OOCC#C(23830)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/OneDe] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cs_H_out_H/Ct]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]#C[CH]OOC=C(23831)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.48073e+07,'s^-1'), n=1.57892, Ea=(169.446,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Y_rad_out;Cd_H_out_doubleC] for rate rule [R6HJ_2;Ct_rad_out;Cd_H_out_doubleC]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['C=C1C=[C][CH]OO1(23832)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.23539e+09,'s^-1'), n=0.445806, Ea=(31.0324,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C=COO[C]=C(22461)'],
    products = ['C#CC1OOC1=C(23833)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]OO[C]=C(4541)', '[C]#C(5143)'],
    products = ['[CH]=C=COO[C]=C(22461)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

network(
    label = '4758',
    isomers = [
        '[CH]=C=COO[C]=C(22461)',
    ],
    reactants = [
        ('C=C=O(598)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4758',
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

