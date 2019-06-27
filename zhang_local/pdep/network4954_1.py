species(
    label = '[O]C=C=CC=C=C[O](22651)',
    structure = SMILES('[O]C=C=CC=C=C[O]'),
    E0 = (240.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,540,563.333,586.667,610,1970,2140,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.4721,'amu*angstrom^2'), symmetry=1, barrier=(33.8464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19862,0.0651468,-2.4668e-05,-3.814e-08,2.60145e-11,29069.7,25.2053], Tmin=(100,'K'), Tmax=(917.411,'K')), NASAPolynomial(coeffs=[26.525,-0.000893243,3.60937e-06,-7.7177e-10,4.87264e-14,22188,-110.715], Tmin=(917.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[O]C=[C]C1C=C=CO1(25116)',
    structure = SMILES('[O]C=[C]C1C=C=CO1'),
    E0 = (524.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731237,0.0470005,2.80629e-05,-9.16092e-08,4.41759e-11,63245.1,23.1692], Tmin=(100,'K'), Tmax=(931.324,'K')), NASAPolynomial(coeffs=[26.9262,-0.00172505,3.8147e-06,-7.17773e-10,3.83377e-14,55599.9,-116.193], Tmin=(931.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1,2-Cyclopentadiene) + radical(Cds_S) + radical(C=COJ)"""),
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
    label = '[O]C=C=CC=C=C=O(25117)',
    structure = SMILES('[O]C=C=CC=C=C=O'),
    E0 = (284.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,563.333,586.667,610,1970,2140,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.266264,'amu*angstrom^2'), symmetry=1, barrier=(6.12192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47826,0.0415898,-6.4162e-06,-3.59309e-08,2.07806e-11,34322.9,8.50683], Tmin=(100,'K'), Tmax=(928.191,'K')), NASAPolynomial(coeffs=[19.2881,0.000626489,1.94831e-06,-3.99681e-10,2.24634e-14,29475.1,-84.3927], Tmin=(928.191,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[O]C=[C]C=CC#CO(25118)',
    structure = SMILES('[O]C=[C]C=CC#CO'),
    E0 = (312.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3615,1277.5,1000,2100,2250,500,550,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.52429,'amu*angstrom^2'), symmetry=1, barrier=(35.0464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52431,'amu*angstrom^2'), symmetry=1, barrier=(35.0468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52448,'amu*angstrom^2'), symmetry=1, barrier=(35.0508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.4566,0.0760461,-7.84218e-05,3.7756e-08,-6.73433e-12,37785.9,29.7427], Tmin=(100,'K'), Tmax=(1590.47,'K')), NASAPolynomial(coeffs=[22.546,0.00550999,6.56237e-08,-1.57737e-10,1.34234e-14,32073.4,-86.8613], Tmin=(1590.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-CtH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C=[C]C=C=C=CO(25119)',
    structure = SMILES('[O]C=[C]C=C=C=CO'),
    E0 = (297.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49843,'amu*angstrom^2'), symmetry=1, barrier=(34.452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49556,'amu*angstrom^2'), symmetry=1, barrier=(34.3858,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59623,0.0940222,-0.000110897,5.85044e-08,-1.11014e-11,36061.8,30.5288], Tmin=(100,'K'), Tmax=(1554.76,'K')), NASAPolynomial(coeffs=[27.0716,-0.00338443,5.89737e-06,-1.36062e-09,9.79775e-14,30006.1,-111.235], Tmin=(1554.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C=C=C=C[C]=CO(25120)',
    structure = SMILES('[O]C=C=C=C[C]=CO'),
    E0 = (297.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,540,563.333,586.667,610,1970,2140,3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.49843,'amu*angstrom^2'), symmetry=1, barrier=(34.452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49556,'amu*angstrom^2'), symmetry=1, barrier=(34.3858,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.59623,0.0940222,-0.000110897,5.85044e-08,-1.11014e-11,36061.8,30.5288], Tmin=(100,'K'), Tmax=(1554.76,'K')), NASAPolynomial(coeffs=[27.0716,-0.00338443,5.89737e-06,-1.36062e-09,9.79775e-14,30006.1,-111.235], Tmin=(1554.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=COJ)"""),
)

species(
    label = 'O=C=[C][CH]C=C=CO(25121)',
    structure = SMILES('O=C=[C][CH]C=C=CO'),
    E0 = (237.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.620168,0.0821118,-9.08641e-05,4.66535e-08,-8.85776e-12,28743.3,27.6213], Tmin=(100,'K'), Tmax=(1484.59,'K')), NASAPolynomial(coeffs=[23.7047,0.0038827,9.98464e-07,-3.5564e-10,2.79892e-14,22919.3,-94.6339], Tmin=(1484.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJC=C=O) + radical(CCCJ=C=O)"""),
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
    label = '[O]C=[C]C=C=C=C[O](25122)',
    structure = SMILES('[O]C=[C]C=C=C=C[O]'),
    E0 = (439.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,563.333,586.667,610,1970,2140,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63651,'amu*angstrom^2'), symmetry=1, barrier=(37.6267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.929533,0.0840019,-9.66015e-05,5.02119e-08,-9.46959e-12,53047.8,29.5487], Tmin=(100,'K'), Tmax=(1544.09,'K')), NASAPolynomial(coeffs=[24.6852,-0.000873238,3.84058e-06,-9.2181e-10,6.69419e-14,47345.2,-98.0068], Tmin=(1544.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C=C=C[CH][C]=C=O(25123)',
    structure = SMILES('[O]C=C=C[CH][C]=C=O'),
    E0 = (378.926,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.69575,'amu*angstrom^2'), symmetry=1, barrier=(38.9887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.69823,'amu*angstrom^2'), symmetry=1, barrier=(39.0456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0479189,0.0720799,-7.65467e-05,3.83537e-08,-7.22883e-12,45729.3,26.6359], Tmin=(100,'K'), Tmax=(1439.8,'K')), NASAPolynomial(coeffs=[21.0316,0.00680077,-1.26314e-06,1.26794e-10,-6.38645e-15,40410.6,-79.7331], Tmin=(1439.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCCJ=C=O) + radical(C=COJ) + radical(C=CCJC=C=O)"""),
)

species(
    label = '[O]C=C=CC=C1[CH]O1(25124)',
    structure = SMILES('[O]C=C=CC=C1[CH]O1'),
    E0 = (266.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.578408,0.0444297,5.09748e-05,-1.28698e-07,6.12376e-11,32216.6,23.117], Tmin=(100,'K'), Tmax=(914.393,'K')), NASAPolynomial(coeffs=[31.8392,-0.0114088,9.84413e-06,-1.93937e-09,1.23403e-13,23117.2,-143.397], Tmin=(914.393,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(methyleneoxirane) + radical(C=COJ) + radical(C=CCJO)"""),
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
    label = '[O]C=C1[CH]C=C=CO1(25126)',
    structure = SMILES('[O]C=C1[CH]C=C=CO1'),
    E0 = (159.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35549,0.0287868,8.18253e-05,-1.38269e-07,5.7176e-11,19302.8,18.3188], Tmin=(100,'K'), Tmax=(957.746,'K')), NASAPolynomial(coeffs=[23.8802,0.00851092,-1.99972e-06,5.32619e-10,-5.61201e-14,11603.6,-107.047], Tmin=(957.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(C=CCJCO) + radical(C=COJ)"""),
)

species(
    label = '[O]C1[C]=CC=C=CO1(25127)',
    structure = SMILES('[O]C1[C]=CC=C=CO1'),
    E0 = (412.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1103,0.049852,-9.93375e-06,-2.87838e-08,1.58135e-11,49755.8,17.7417], Tmin=(100,'K'), Tmax=(985.188,'K')), NASAPolynomial(coeffs=[17.4475,0.0153057,-5.72954e-06,1.11913e-09,-8.46981e-14,44994.2,-68.6604], Tmin=(985.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = 'O=C=C=CC=C=CO(25128)',
    structure = SMILES('O=C=C=CC=C=CO'),
    E0 = (143.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06126,0.0486105,-9.91014e-06,-4.21982e-08,2.56428e-11,17326.1,8.59426], Tmin=(100,'K'), Tmax=(910.111,'K')), NASAPolynomial(coeffs=[22.6336,-0.00321583,4.66063e-06,-9.75322e-10,6.37855e-14,11619.1,-103.229], Tmin=(910.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds)"""),
)

species(
    label = '[O]C=[C]C[C]=C=C[O](25129)',
    structure = SMILES('[O]C=[C]C[C]=C=C[O]'),
    E0 = (570.352,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21092,'amu*angstrom^2'), symmetry=1, barrier=(27.8413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21173,'amu*angstrom^2'), symmetry=1, barrier=(27.86,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.348659,0.0777277,-8.37973e-05,4.22072e-08,-7.91547e-12,68769.4,32.7862], Tmin=(100,'K'), Tmax=(1482.05,'K')), NASAPolynomial(coeffs=[22.7391,0.00518801,-2.85037e-08,-1.30243e-10,1.16247e-14,63049.1,-83.8929], Tmin=(1482.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]C=[C]C=C[O](25130)',
    structure = SMILES('[O]C=[C]C=[C]C=C[O]'),
    E0 = (409.024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63376,'amu*angstrom^2'), symmetry=1, barrier=(37.5634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62379,'amu*angstrom^2'), symmetry=1, barrier=(37.3341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.512,0.0895128,-0.00010131,5.21083e-08,-9.62984e-12,49420.1,32.1572], Tmin=(100,'K'), Tmax=(1614.26,'K')), NASAPolynomial(coeffs=[24.8535,-0.000370697,5.02587e-06,-1.22875e-09,8.96403e-14,44106.8,-97.7586], Tmin=(1614.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C=[C][CH]C=C[C]=O(25131)',
    structure = SMILES('[O]C=[C][CH]C=C[C]=O'),
    E0 = (398.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.2047,'amu*angstrom^2'), symmetry=1, barrier=(27.6984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20479,'amu*angstrom^2'), symmetry=1, barrier=(27.7005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20531,'amu*angstrom^2'), symmetry=1, barrier=(27.7125,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.455792,0.0687658,-6.37779e-05,2.49092e-08,-2.59275e-12,48018.3,29.5609], Tmin=(100,'K'), Tmax=(1046.04,'K')), NASAPolynomial(coeffs=[18.7572,0.0133945,-5.33013e-06,1.01314e-09,-7.32685e-14,43390.1,-63.3777], Tmin=(1046.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=COJ) + radical(C=CCJC=C) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C][C]=CC=C[O](25132)',
    structure = SMILES('[O]C=[C][C]=CC=C[O]'),
    E0 = (409.024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63376,'amu*angstrom^2'), symmetry=1, barrier=(37.5634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62379,'amu*angstrom^2'), symmetry=1, barrier=(37.3341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.512,0.0895128,-0.00010131,5.21083e-08,-9.62984e-12,49420.1,32.1572], Tmin=(100,'K'), Tmax=(1614.26,'K')), NASAPolynomial(coeffs=[24.8535,-0.000370697,5.02587e-06,-1.22875e-09,8.96403e-14,44106.8,-97.7586], Tmin=(1614.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJC=C) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C[C][CH]C#CC=O(25133)',
    structure = SMILES('[O]C[C][CH]C#CC=O'),
    E0 = (656.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,2100,2250,500,550,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863454,0.0773707,-0.000123236,1.09372e-07,-3.75205e-11,79020.8,29.5587], Tmin=(100,'K'), Tmax=(875.041,'K')), NASAPolynomial(coeffs=[6.65758,0.0324874,-1.47604e-05,2.70184e-09,-1.80235e-13,78711.1,6.40333], Tmin=(875.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCJ2_triplet) + radical(CCOJ) + radical(Sec_Propargyl)"""),
)

species(
    label = '[O]C=[C]C[CH][C]=C=O(25134)',
    structure = SMILES('[O]C=[C]C[CH][C]=C=O'),
    E0 = (498.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,340.729,340.839,341.32,342.49],'cm^-1')),
        HinderedRotor(inertia=(0.31352,'amu*angstrom^2'), symmetry=1, barrier=(25.8298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310355,'amu*angstrom^2'), symmetry=1, barrier=(25.8262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312533,'amu*angstrom^2'), symmetry=1, barrier=(25.8327,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646728,0.0717947,-7.92167e-05,4.40599e-08,-9.65458e-12,60094.2,30.0041], Tmin=(100,'K'), Tmax=(1113.39,'K')), NASAPolynomial(coeffs=[15.3167,0.0190908,-8.21167e-06,1.54384e-09,-1.08003e-13,56827.5,-42.3455], Tmin=(1113.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(CCCJ=C=O) + radical(Cds_S) + radical(CCJC(C)=C=O) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C[CH][CH][C]=C=O(25135)',
    structure = SMILES('[O]C=C[CH][CH][C]=C=O'),
    E0 = (401.904,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,437.754,437.87,438.12,438.583],'cm^-1')),
        HinderedRotor(inertia=(0.257009,'amu*angstrom^2'), symmetry=1, barrier=(35.0878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256834,'amu*angstrom^2'), symmetry=1, barrier=(35.1163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.256811,'amu*angstrom^2'), symmetry=1, barrier=(35.1012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660128,0.0645532,-5.45706e-05,1.81813e-08,-8.18405e-13,48465.9,28.3742], Tmin=(100,'K'), Tmax=(1039.48,'K')), NASAPolynomial(coeffs=[17.1877,0.0160285,-6.30108e-06,1.17505e-09,-8.36442e-14,44215.4,-55.9191], Tmin=(1039.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(Allyl_S) + radical(CCJC(C)=C=O) + radical(C=COJ) + radical(CCCJ=C=O)"""),
)

species(
    label = '[O]C[C]=C[CH][C]=C=O(25136)',
    structure = SMILES('[O]C[C]=C[CH][C]=C=O'),
    E0 = (580.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,206.104,206.104,206.104,206.104],'cm^-1')),
        HinderedRotor(inertia=(0.0304393,'amu*angstrom^2'), symmetry=1, barrier=(55.0712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00663462,'amu*angstrom^2'), symmetry=1, barrier=(12.0034,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.82695,'amu*angstrom^2'), symmetry=1, barrier=(55.0712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00165,0.0732859,-0.000111412,9.92961e-08,-3.53494e-11,69866.4,27.9704], Tmin=(100,'K'), Tmax=(819.267,'K')), NASAPolynomial(coeffs=[6.3683,0.0343176,-1.66906e-05,3.19821e-09,-2.21028e-13,69415.4,5.76394], Tmin=(819.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + radical(Cds_S) + radical(C=CCJC=C=O) + radical(CCCJ=C=O) + radical(CCOJ)"""),
)

species(
    label = 'O=C=[C][CH]C=CC=O(25137)',
    structure = SMILES('O=C=[C][CH]C=CC=O'),
    E0 = (182.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32895,0.0627944,-6.0756e-05,3.06259e-08,-6.42512e-12,22005.1,22.9869], Tmin=(100,'K'), Tmax=(1115.78,'K')), NASAPolynomial(coeffs=[10.9653,0.0282483,-1.43132e-05,2.87645e-09,-2.07543e-13,19854.8,-24.5584], Tmin=(1115.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJC=C=O) + radical(CCCJ=C=O)"""),
)

species(
    label = 'C#CC([O])C([O])C#C(22644)',
    structure = SMILES('C#CC([O])C([O])C#C'),
    E0 = (510.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,756.667,763.333,770,3350,3450,2000,2200,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2100,2250,500,550,238.442,238.442,238.443],'cm^-1')),
        HinderedRotor(inertia=(0.918725,'amu*angstrom^2'), symmetry=1, barrier=(37.0661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.918728,'amu*angstrom^2'), symmetry=1, barrier=(37.0661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.918718,'amu*angstrom^2'), symmetry=1, barrier=(37.0661,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0752804,0.0763947,-8.55937e-05,4.621e-08,-9.34535e-12,61582.1,30.3555], Tmin=(100,'K'), Tmax=(1369.22,'K')), NASAPolynomial(coeffs=[20.1033,0.00816174,-6.73045e-07,-8.94211e-11,1.24199e-14,56926.6,-70.1578], Tmin=(1369.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Ct-CtCs) + group(Ct-CtCs) + group(Ct-CtH) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[O]C=C1C=C[C]1C=O(25138)',
    structure = SMILES('[O]C=C1C=C[C]1C=O'),
    E0 = (127.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29941,0.0357901,4.84943e-05,-1.03081e-07,4.59264e-11,15506.6,19.5471], Tmin=(100,'K'), Tmax=(934.994,'K')), NASAPolynomial(coeffs=[22.454,0.00545408,6.3934e-07,-1.36787e-10,-7.05683e-16,8920.83,-95.1533], Tmin=(934.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsCsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(C=COJ) + radical(C=CCJ(C=O)C=C)"""),
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
    label = '[O]C=[C]C1C=C1C=O(25077)',
    structure = SMILES('[O]C=[C]C1C=C1C=O'),
    E0 = (390.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928459,0.0665916,-6.76058e-05,3.42496e-08,-6.91644e-12,47085.3,24.2779], Tmin=(100,'K'), Tmax=(1191.7,'K')), NASAPolynomial(coeffs=[14.5976,0.020711,-9.85651e-06,1.94368e-09,-1.3926e-13,43827.3,-44.0655], Tmin=(1191.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopropene) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CC#CC=O(25139)',
    structure = SMILES('[O]C=C=CC#CC=O'),
    E0 = (213.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2100,2250,500,550,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.52758,'amu*angstrom^2'), symmetry=1, barrier=(35.122,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52746,'amu*angstrom^2'), symmetry=1, barrier=(35.1192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.801567,0.0595502,-4.30813e-05,3.96712e-09,4.86838e-12,25862.4,25.326], Tmin=(100,'K'), Tmax=(998.145,'K')), NASAPolynomial(coeffs=[18.6072,0.010938,-4.20446e-06,8.2813e-10,-6.27893e-14,21175,-66.2177], Tmin=(998.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=C[C]=CC=O(25140)',
    structure = SMILES('[O]C=C=C[C]=CC=O'),
    E0 = (242.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.34869,0.0747965,-8.12054e-05,4.29916e-08,-8.8381e-12,29323.5,25.9053], Tmin=(100,'K'), Tmax=(1193.62,'K')), NASAPolynomial(coeffs=[18.1178,0.0152491,-6.37297e-06,1.19558e-09,-8.40089e-14,25081.6,-62.9655], Tmin=(1193.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]C=C=C[CH]C=C=O(25141)',
    structure = SMILES('[O]C=C=C[CH]C=C=O'),
    E0 = (176.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.559273,0.0653097,-4.98482e-05,7.07922e-09,4.87741e-12,21378.8,25.0833], Tmin=(100,'K'), Tmax=(964.311,'K')), NASAPolynomial(coeffs=[19.0752,0.0121245,-3.85842e-06,6.84982e-10,-4.99716e-14,16709.6,-69.2668], Tmin=(964.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CCJC=C=O)"""),
)

species(
    label = '[O]C=C=C=CC=C[O](25142)',
    structure = SMILES('[O]C=C=C=CC=C[O]'),
    E0 = (240.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.19862,0.0651468,-2.4668e-05,-3.814e-08,2.60145e-11,29069.7,25.8985], Tmin=(100,'K'), Tmax=(917.411,'K')), NASAPolynomial(coeffs=[26.525,-0.000893243,3.60937e-06,-7.7177e-10,4.87264e-14,22188,-110.021], Tmin=(917.411,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C1[C]=CC=C1C=O(25143)',
    structure = SMILES('[O]C1[C]=CC=C1C=O'),
    E0 = (314.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30802,0.0485238,-1.99353e-05,-8.21868e-09,5.68045e-12,37973.2,24.1116], Tmin=(100,'K'), Tmax=(1131.26,'K')), NASAPolynomial(coeffs=[15.1175,0.0203413,-9.94282e-06,2.0258e-09,-1.48833e-13,33527.7,-50.0535], Tmin=(1131.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclopentadiene) + radical(CC(C)OJ) + radical(1,3-cyclopentadiene-vinyl-1)"""),
)

species(
    label = '[C]1=COOC=[C]C=C1(25144)',
    structure = SMILES('[C]1=COOC=[C]C=C1'),
    E0 = (569.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72199,0.0246986,7.91798e-05,-1.34129e-07,5.73625e-11,68558,18.041], Tmin=(100,'K'), Tmax=(920.683,'K')), NASAPolynomial(coeffs=[20.9522,0.00678973,1.41697e-06,-3.85229e-10,1.94154e-14,62235.1,-88.2525], Tmin=(920.683,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.145,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclooctane) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = 'O=C=C=CC=CC=O(25145)',
    structure = SMILES('O=C=C=CC=CC=O'),
    E0 = (87.7784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58428,0.0458329,-3.64185e-05,1.3618e-08,-1.98663e-12,10650.3,9.0951], Tmin=(100,'K'), Tmax=(1640,'K')), NASAPolynomial(coeffs=[15.4047,0.0121242,-5.58697e-06,1.0847e-09,-7.60309e-14,6117.31,-64.4169], Tmin=(1640,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.7784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds)"""),
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
    label = '[C]=CC=C=C[O](23395)',
    structure = SMILES('[C]=CC=C=C[O]'),
    E0 = (725.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.41438,'amu*angstrom^2'), symmetry=1, barrier=(32.5194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.675917,0.0557995,-5.93347e-05,2.89865e-08,-5.19467e-12,87364.6,22.4469], Tmin=(100,'K'), Tmax=(1604.96,'K')), NASAPolynomial(coeffs=[17.9233,0.00199229,1.06811e-06,-3.04893e-10,2.23918e-14,83222.2,-64.5789], Tmin=(1604.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(725.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'O=CC1=CC=C1C=O(25074)',
    structure = SMILES('O=CC1=CC=C1C=O'),
    E0 = (195.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.559466,0.0640819,-5.04012e-05,1.32407e-08,3.19995e-13,23661,20.9059], Tmin=(100,'K'), Tmax=(1134.84,'K')), NASAPolynomial(coeffs=[19.873,0.0136584,-7.08416e-06,1.50005e-09,-1.1301e-13,18140.8,-79.7216], Tmin=(1134.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + ring(cyclobutadiene_13)"""),
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
    E0 = (240.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (524.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (499.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (490.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (592.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (330.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (312.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (539.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (651.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (590.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (427.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (373.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (303.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (412.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (265.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (598.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (431.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (425.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (497.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (695.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (516.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (426.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (605.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (404.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (655.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (362.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (861.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (390.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (439.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (396.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (427.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (403.198,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (462.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (315.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (569.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (279.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (792.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (248.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['C#CC=O(21959)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=[C]C1C=C=CO1(25116)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2e+11,'s^-1'), n=0.21, Ea=(284.262,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMM;doublebond_intra;radadd_intra] for rate rule [R6_SMM;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 281.8 to 284.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[O]C=C=CC=C=C=O(25117)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C=[C]C=CC#CO(25118)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C=[C]C=C=C=CO(25119)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(8.63625e+10,'s^-1'), n=1.0925, Ea=(294.328,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R4H_MMS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C=C=C=C[C]=CO(25120)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out;XH_out] for rate rule [R5H_SMMS;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.73205080757
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['O=C=[C][CH]C=C=CO(25121)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.25419e+06,'s^-1'), n=1.03067, Ea=(72.2138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7H;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C=C[O](8556)', '[CH]=C=C[O](8556)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.145e+09,'cm^3/(mol*s)'), n=0.8, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 131 used for Cd_allenic;Cd_allenic
Exact match found for rate rule [Cd_allenic;Cd_allenic]
Euclidian distance = 0
family: R_Recombination
Ea raised from -4.3 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(8)', '[O]C=[C]C=C=C=C[O](25122)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(8)', '[O]C=C=C[CH][C]=C=O(25123)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=C=CC=C1[CH]O1(25124)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.31909e+12,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=C=CC1[C]=CO1(25125)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.94552e+11,'s^-1'), n=0.183155, Ea=(133.42,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 131.2 to 133.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=C1[CH]C=C=CO1(25126)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.9903e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C1[C]=CC=C=CO1(25127)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.44e+12,'s^-1'), n=-0.622, Ea=(172.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;doublebond_intra_CdCdd;radadd_intra] for rate rule [R7;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 169.8 to 172.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['O=C=C=CC=C=CO(25128)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C=[C]C[C]=C=C[O](25129)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C=[C]C=[C]C=C[O](25130)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O]C=[C][CH]C=C[C]=O(25131)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O]C=[C][C]=CC=C[O](25132)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C[C][CH]C#CC=O(25133)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C=[C]C[CH][C]=C=O(25134)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C=C[CH][CH][C]=C=O(25135)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C[C]=C[CH][C]=C=O(25136)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['O=C=[C][CH]C=CC=O(25137)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.16796e+09,'s^-1'), n=0.809263, Ea=(163.807,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [1_3_4_pentatriene;CH_end;CddC_2]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_ene_reaction"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C#CC([O])C([O])C#C(22644)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.30946e+10,'s^-1'), n=0.360276, Ea=(144.706,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for 1_5_hexadiyne
Exact match found for rate rule [1_5_hexadiyne]
Euclidian distance = 0
family: 6_membered_central_C-C_shift"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=C1C=C[C]1C=O(25138)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.99998e+11,'s^-1'), n=0.0559095, Ea=(122.413,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for 1,3-butadiene_backbone;CddC_1;CddC_2
Exact match found for rate rule [1,3-butadiene_backbone;CddC_1;CddC_2]
Euclidian distance = 0
family: Intra_2+2_cycloaddition_Cd"""),
)

reaction(
    label = 'reaction27',
    reactants = ['O(T)(63)', 'C#CC=C[C]=C[O](22635)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=[C]C1C=C1C=O(25077)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.252e+14,'s^-1'), n=-0.355, Ea=(150.144,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 150.0 to 150.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(8)', '[O]C=C=CC#CC=O(25139)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C=C[O](8556)', 'C#CC=O(21959)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.5e-18,'cm^3/(molecule*s)'), n=1.925, Ea=(42.091,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CdsJ=Cdd] for rate rule [Ct-H_Ct-CO;CdsJ=Cdd]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=C=C[C]=CC=O(25140)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.55304e+08,'s^-1'), n=1.65613, Ea=(187.464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_D;Cd_rad_out_single;Cd_H_out_singleDe] + [R2H_D;Cd_rad_out_singleDe;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe;Cd_H_out_singleDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=C=C[CH]C=C=O(25141)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(9.93038e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C=C=C=CC=C[O](25142)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(383,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_single;Cd_H_out_doubleC] for rate rule [R3H_DS;Cd_rad_out_singleDe;Cd_H_out_doubleC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[O]C1[C]=CC=C1C=O(25143)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(9.54e+12,'s^-1'), n=-0.453, Ea=(75.061,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 833 used for R5_MS;doublebond_intra_CdCdd;radadd_intra_cdsingleDe
Exact match found for rate rule [R5_MS;doublebond_intra_CdCdd;radadd_intra_cdsingleDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['[C]1=COOC=[C]C=C1(25144)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(328.727,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra] for rate rule [R8_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 323.3 to 328.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['O=C=C=CC=CC=O(25145)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.70659e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=O(373)', '[C]=CC=C=C[O](23395)'],
    products = ['[O]C=C=CC=C=C[O](22651)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C=C=CC=C=C[O](22651)'],
    products = ['O=CC1=CC=C1C=O(25074)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleDe_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

network(
    label = '4954',
    isomers = [
        '[O]C=C=CC=C=C[O](22651)',
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
    label = '4954',
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

