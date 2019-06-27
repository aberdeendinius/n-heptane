species(
    label = '[O]C=C=CC([O])O[O](22518)',
    structure = SMILES('[O]C=C=CC([O])O[O]'),
    E0 = (168.529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.790077,'amu*angstrom^2'), symmetry=1, barrier=(18.1654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.790635,'amu*angstrom^2'), symmetry=1, barrier=(18.1782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.820596,0.075151,-0.00010872,8.4201e-08,-2.61703e-11,20379,28.9596], Tmin=(100,'K'), Tmax=(786.749,'K')), NASAPolynomial(coeffs=[10.7806,0.0245088,-1.216e-05,2.37367e-09,-1.66865e-13,18811.9,-16.7021], Tmin=(786.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC=O(5472)',
    structure = SMILES('[O]OC=O'),
    E0 = (-195.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,739.225,739.248,739.254,739.261,739.261],'cm^-1')),
        HinderedRotor(inertia=(0.00030847,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3570.08,'J/mol'), sigma=(5.61676,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=557.64 K, Pc=45.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.38375,0.00191416,4.59985e-05,-6.61808e-08,2.67351e-11,-23479.8,12.2589], Tmin=(100,'K'), Tmax=(935.456,'K')), NASAPolynomial(coeffs=[10.7708,0.000103189,1.15693e-06,-1.97269e-10,7.46894e-15,-26164.6,-29.8498], Tmin=(935.456,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cds-OdOsH) + radical(C(=O)OOJ)"""),
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
    label = '[O]C=[C]C1OC1O[O](24581)',
    structure = SMILES('[O]C=[C]C1OC1O[O]'),
    E0 = (194.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.162914,0.0716953,-7.96721e-05,4.2252e-08,-8.16187e-12,23502.9,31.0528], Tmin=(100,'K'), Tmax=(1528.78,'K')), NASAPolynomial(coeffs=[18.6555,0.00540026,2.11128e-06,-7.10329e-10,5.72177e-14,19742.3,-61.2037], Tmin=(1528.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(ROOJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]C1OOC1[O](24670)',
    structure = SMILES('[O]C=[C]C1OOC1[O]'),
    E0 = (253.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0344073,0.0644751,-6.31693e-05,3.00699e-08,-5.25491e-12,30652.7,31.8249], Tmin=(100,'K'), Tmax=(1695.56,'K')), NASAPolynomial(coeffs=[16.9965,0.00860276,1.43505e-07,-2.82813e-10,2.53518e-14,27133.3,-52.6794], Tmin=(1695.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(253.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(12dioxetane) + radical(CCOJ) + radical(Cds_S) + radical(C=COJ)"""),
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
    label = '[O]C=C=CC(=O)O[O](24671)',
    structure = SMILES('[O]C=C=CC(=O)O[O]'),
    E0 = (11.6317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,338.88,338.88,338.881,338.884,338.894,338.896],'cm^-1')),
        HinderedRotor(inertia=(0.215091,'amu*angstrom^2'), symmetry=1, barrier=(17.5279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.21509,'amu*angstrom^2'), symmetry=1, barrier=(17.5278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.715101,0.065582,-7.5662e-05,4.15083e-08,-8.69759e-12,1522.96,28.8817], Tmin=(100,'K'), Tmax=(1183.32,'K')), NASAPolynomial(coeffs=[17.6998,0.00816784,-2.88221e-06,5.04826e-10,-3.47327e-14,-2496.68,-55.9185], Tmin=(1183.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.6317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C(=O)OOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]OC([O])C=C=C=O(24672)',
    structure = SMILES('[O]OC([O])C=C=C=O'),
    E0 = (212.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,540,610,2055,1380,1390,370,380,2900,435,2120,512.5,787.5,3010,987.5,1337.5,450,1655,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000570816,'amu*angstrom^2'), symmetry=1, barrier=(6.48107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000570135,'amu*angstrom^2'), symmetry=1, barrier=(6.47335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9564,0.0534695,-9.79374e-05,9.75516e-08,-3.68828e-11,25638.2,12.0727], Tmin=(100,'K'), Tmax=(831.447,'K')), NASAPolynomial(coeffs=[3.88852,0.025405,-1.34454e-05,2.65419e-09,-1.85353e-13,25965.7,7.00939], Tmin=(831.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCOJ)"""),
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
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62683,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1487.4],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,-1038.59,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,-1040.82,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62683,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[O]C=[C]C=C(O)O[O](24673)',
    structure = SMILES('[O]C=[C]C=C(O)O[O]'),
    E0 = (136.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.196147,0.0881846,-0.000124526,8.22469e-08,-2.03272e-11,16530.6,28.3643], Tmin=(100,'K'), Tmax=(1092.77,'K')), NASAPolynomial(coeffs=[20.3419,0.00524725,-3.01484e-08,-2.02226e-10,2.18089e-14,12505.2,-70.4219], Tmin=(1092.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O]OC([O])C=C=[C]O(24674)',
    structure = SMILES('[O]OC([O])C=C=[C]O'),
    E0 = (266.81,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,540,610,2055,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.569363,'amu*angstrom^2'), symmetry=1, barrier=(13.0908,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569014,'amu*angstrom^2'), symmetry=1, barrier=(13.0828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.568475,'amu*angstrom^2'), symmetry=1, barrier=(13.0704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.46464,0.0858319,-0.000146463,1.27323e-07,-4.25345e-11,32209.5,31.9467], Tmin=(100,'K'), Tmax=(859.545,'K')), NASAPolynomial(coeffs=[10.5063,0.0238682,-1.17448e-05,2.21691e-09,-1.5006e-13,31045.9,-11.7052], Tmin=(859.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=CJO) + radical(ROOJ)"""),
)

species(
    label = '[O]C=[C][CH]C(=O)OO(24675)',
    structure = SMILES('[O]C=[C][CH]C(=O)OO'),
    E0 = (-73.5576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.511965,0.0571377,-1.05058e-05,-4.0751e-08,2.21476e-11,-8703.89,26.6669], Tmin=(100,'K'), Tmax=(990.642,'K')), NASAPolynomial(coeffs=[24.3127,0.00666199,-3.16327e-06,8.00249e-10,-7.12172e-14,-15658.3,-99.2339], Tmin=(990.642,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-73.5576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=[C]C(O)O[O](24676)',
    structure = SMILES('[O]C=C=[C]C(O)O[O]'),
    E0 = (180.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,540,610,2055,1380,1390,370,380,2900,435,1685,370,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.783949,'amu*angstrom^2'), symmetry=1, barrier=(18.0245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.785787,'amu*angstrom^2'), symmetry=1, barrier=(18.0668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.78308,'amu*angstrom^2'), symmetry=1, barrier=(18.0046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.395181,0.082893,-0.000122773,9.06153e-08,-2.60486e-11,21855.6,30.097], Tmin=(100,'K'), Tmax=(858.573,'K')), NASAPolynomial(coeffs=[14.4921,0.0172159,-8.02747e-06,1.51656e-09,-1.04419e-13,19434.9,-35.7629], Tmin=(858.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]C=C=[C]C([O])OO(24677)',
    structure = SMILES('[O]C=C=[C]C([O])OO'),
    E0 = (254.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646174,0.0815067,-0.000120009,8.78923e-08,-2.3253e-11,30706.7,29.3224], Tmin=(100,'K'), Tmax=(635.623,'K')), NASAPolynomial(coeffs=[11.2407,0.0255795,-1.33834e-05,2.65389e-09,-1.87609e-13,29142.8,-18.6966], Tmin=(635.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_S) + radical(C=COJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])C#C[CH]O(24678)',
    structure = SMILES('[O]OC([O])C#C[CH]O'),
    E0 = (246.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,3615,1277.5,1000,2100,2250,500,550,1380,1390,370,380,2900,435,180,1354.68,1357.46],'cm^-1')),
        HinderedRotor(inertia=(0.349286,'amu*angstrom^2'), symmetry=1, barrier=(8.03077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.350909,'amu*angstrom^2'), symmetry=1, barrier=(8.06809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.351318,'amu*angstrom^2'), symmetry=1, barrier=(8.0775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.6323,'amu*angstrom^2'), symmetry=1, barrier=(60.5218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654606,0.0872388,-0.000164914,1.54312e-07,-5.30816e-11,29720,31.7589], Tmin=(100,'K'), Tmax=(911.509,'K')), NASAPolynomial(coeffs=[6.13414,0.0291168,-1.31905e-05,2.3295e-09,-1.48556e-13,30136.7,13.5962], Tmin=(911.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(CCsJOH) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC(O)[CH][C]=C=O(24679)',
    structure = SMILES('[O]OC(O)[CH][C]=C=O'),
    E0 = (80.9653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456576,0.0834113,-0.000126337,9.93321e-08,-3.09561e-11,9860.38,27.8305], Tmin=(100,'K'), Tmax=(787.915,'K')), NASAPolynomial(coeffs=[12.3324,0.0231289,-1.15883e-05,2.25319e-09,-1.57465e-13,7988.72,-26.6342], Tmin=(787.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.9653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + radical(ROOJ) + radical(C=CCJCO) + radical(CCCJ=C=O)"""),
)

species(
    label = '[O]OC(=O)[CH][C]=CO(24680)',
    structure = SMILES('[O]OC(=O)[CH][C]=CO'),
    E0 = (-20.6316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.377932,0.0603949,-1.83972e-05,-3.75969e-08,2.29356e-11,-2333.51,26.4946], Tmin=(100,'K'), Tmax=(960.699,'K')), NASAPolynomial(coeffs=[25.6238,0.00191147,1.08746e-07,8.5371e-11,-1.81608e-14,-9336.13,-105.489], Tmin=(960.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.6316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(Cds_S) + radical(C=CCJCO) + radical(C(=O)OOJ)"""),
)

species(
    label = '[O]C([CH][C]=C=O)OO(24681)',
    structure = SMILES('[O]C([CH][C]=C=O)OO'),
    E0 = (154.666,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.484043,0.0852237,-0.000138064,1.21895e-07,-4.30217e-11,18721,27.8256], Tmin=(100,'K'), Tmax=(777.188,'K')), NASAPolynomial(coeffs=[9.31565,0.0310473,-1.66671e-05,3.32142e-09,-2.34688e-13,17611.7,-10.8608], Tmin=(777.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCCJ=C=O) + radical(C=CCJCO) + radical(CCOJ)"""),
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
    label = '[O][CH]O[O](8201)',
    structure = SMILES('[O][CH]O[O]'),
    E0 = (228.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,1824.75],'cm^-1')),
        HinderedRotor(inertia=(0.270955,'amu*angstrom^2'), symmetry=1, barrier=(6.2298,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (61.0168,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.97073,0.0345271,-8.63149e-05,9.83474e-08,-3.87219e-11,27554.6,13.3771], Tmin=(100,'K'), Tmax=(878.005,'K')), NASAPolynomial(coeffs=[-0.738633,0.0210949,-1.15487e-05,2.23205e-09,-1.51294e-13,29375,37.4477], Tmin=(878.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-OsOsHH) + radical(OCOJ) + radical(ROOJ) + radical(OCJO)"""),
)

species(
    label = '[O]C=[C][CH]C(=O)O[O](24682)',
    structure = SMILES('[O]C=[C][CH]C(=O)O[O]'),
    E0 = (120.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796009,0.0533235,-1.45159e-05,-3.22026e-08,1.86425e-11,14663.3,26.4055], Tmin=(100,'K'), Tmax=(986.773,'K')), NASAPolynomial(coeffs=[22.4406,0.0054861,-2.45264e-06,6.25951e-10,-5.66116e-14,8448.99,-87.5724], Tmin=(986.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(C(=O)OOJ) + radical(C=CCJCO) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=[C]C([O])O[O](24683)',
    structure = SMILES('[O]C=C=[C]C([O])O[O]'),
    E0 = (406.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.706711,'amu*angstrom^2'), symmetry=1, barrier=(16.2487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706741,'amu*angstrom^2'), symmetry=1, barrier=(16.2494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654898,0.0814635,-0.000140229,1.23307e-07,-4.18647e-11,48988,30.0184], Tmin=(100,'K'), Tmax=(841.746,'K')), NASAPolynomial(coeffs=[10.0037,0.0232706,-1.19961e-05,2.31558e-09,-1.59132e-13,47901.8,-10.5766], Tmin=(841.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCOJ) + radical(ROOJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC([O])[CH][C]=C=O(24684)',
    structure = SMILES('[O]OC([O])[CH][C]=C=O'),
    E0 = (306.67,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,2120,512.5,787.5,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,180,1732.38],'cm^-1')),
        HinderedRotor(inertia=(1.85646,'amu*angstrom^2'), symmetry=1, barrier=(42.6838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85702,'amu*angstrom^2'), symmetry=1, barrier=(42.6966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.85778,'amu*angstrom^2'), symmetry=1, barrier=(42.714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.705271,0.082138,-0.000144472,1.33105e-07,-4.73284e-11,36993.2,27.7897], Tmin=(100,'K'), Tmax=(831.702,'K')), NASAPolynomial(coeffs=[7.89945,0.0290809,-1.5494e-05,3.03669e-09,-2.10848e-13,36434.9,-1.75456], Tmin=(831.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCOJ) + radical(CCCJ=C=O) + radical(C=CCJCO) + radical(ROOJ)"""),
)

species(
    label = '[O]OC([O])C=C1[CH]O1(24685)',
    structure = SMILES('[O]OC([O])C=C1[CH]O1'),
    E0 = (194.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.952218,0.0574747,-4.41778e-05,8.53263e-09,2.56639e-12,23536.4,27.0602], Tmin=(100,'K'), Tmax=(1014.85,'K')), NASAPolynomial(coeffs=[17.2959,0.0119145,-4.71054e-06,9.16316e-10,-6.79748e-14,19248.1,-56.814], Tmin=(1014.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(methyleneoxirane) + radical(ROOJ) + radical(C=CCJO) + radical(CCOJ)"""),
)

species(
    label = '[O]C=C1[CH]C(O[O])O1(24686)',
    structure = SMILES('[O]C=C1[CH]C(O[O])O1'),
    E0 = (17.6279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538696,0.0545163,4.33671e-06,-7.20238e-08,3.98734e-11,2265.14,21.1837], Tmin=(100,'K'), Tmax=(899.743,'K')), NASAPolynomial(coeffs=[27.3992,-0.00607809,7.29523e-06,-1.55746e-09,1.05155e-13,-4949.22,-118.796], Tmin=(899.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.6279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=CCJCO) + radical(ROOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]OC([O])C1[C]=CO1(24687)',
    structure = SMILES('[O]OC([O])C1[C]=CO1'),
    E0 = (289.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799084,0.0590915,-4.15817e-05,1.12858e-09,6.24983e-12,34944.3,27.984], Tmin=(100,'K'), Tmax=(987.431,'K')), NASAPolynomial(coeffs=[19.2121,0.00921717,-3.36237e-06,6.72986e-10,-5.25504e-14,30103.1,-66.7161], Tmin=(987.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[O]C=C1[CH]C([O])OO1(24688)',
    structure = SMILES('[O]C=C1[CH]C([O])OO1'),
    E0 = (54.8226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47209,0.0374682,2.45255e-05,-6.4809e-08,2.84062e-11,6700.9,25.0985], Tmin=(100,'K'), Tmax=(980.028,'K')), NASAPolynomial(coeffs=[18.332,0.0126753,-4.905e-06,1.04491e-09,-8.48005e-14,1282.23,-66.6863], Tmin=(980.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(54.8226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OC1C=[C]C([O])O1(24641)',
    structure = SMILES('[O]OC1C=[C]C([O])O1'),
    E0 = (193.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72137,0.0557374,-6.39757e-05,4.61683e-08,-1.48066e-11,23364.8,24.0727], Tmin=(100,'K'), Tmax=(734.266,'K')), NASAPolynomial(coeffs=[5.9535,0.0326818,-1.6875e-05,3.40275e-09,-2.45534e-13,22743.3,4.96245], Tmin=(734.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(25dihydrofuran) + radical(ROOJ) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[C]=CC([O])OO1(24689)',
    structure = SMILES('[O]C1[C]=CC([O])OO1'),
    E0 = (243.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60471,0.0410576,-2.36633e-05,5.59032e-09,-4.95341e-13,29369.3,19.7867], Tmin=(100,'K'), Tmax=(2607.26,'K')), NASAPolynomial(coeffs=[21.3811,0.0122513,-7.0906e-06,1.35275e-09,-8.90163e-14,19578.3,-88.7919], Tmin=(2607.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro12dioxin) + radical(CCOJ) + radical(CCOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=C=CC(=O)OO(22529)',
    structure = SMILES('[O]C=C=CC(=O)OO'),
    E0 = (-182.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478148,0.0688529,-6.98198e-05,3.07027e-08,-4.28822e-12,-21846.3,28.9735], Tmin=(100,'K'), Tmax=(1038.56,'K')), NASAPolynomial(coeffs=[19.2602,0.00985649,-3.88159e-06,7.46153e-10,-5.48251e-14,-26467.1,-65.8141], Tmin=(1038.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-182.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[O]OC(=O)C=C=CO(24690)',
    structure = SMILES('[O]OC(=O)C=C=CO'),
    E0 = (-129.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.044103,0.0756452,-9.00606e-05,4.98572e-08,-1.03224e-11,-15462.9,29.8776], Tmin=(100,'K'), Tmax=(1301.53,'K')), NASAPolynomial(coeffs=[21.3072,0.00384818,1.1711e-07,-1.40452e-10,1.2472e-14,-20451.5,-76.2097], Tmin=(1301.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C(=O)OOJ)"""),
)

species(
    label = '[O]OC(O)C=C=C=O(24691)',
    structure = SMILES('[O]OC(O)C=C=C=O'),
    E0 = (-13.0799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7771,0.0537187,-7.49903e-05,5.50692e-08,-1.52194e-11,-1497.6,11.876], Tmin=(100,'K'), Tmax=(648.208,'K')), NASAPolynomial(coeffs=[8.27175,0.0195517,-9.6029e-06,1.88676e-09,-1.33376e-13,-2463.75,-17.599], Tmin=(648.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-13.0799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = '[O]C(C=C=C=O)OO(24692)',
    structure = SMILES('[O]C(C=C=C=O)OO'),
    E0 = (60.6206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71985,0.0567597,-9.23723e-05,8.7659e-08,-3.32599e-11,7366.65,12.1622], Tmin=(100,'K'), Tmax=(776.878,'K')), NASAPolynomial(coeffs=[5.32841,0.0273279,-1.45921e-05,2.93244e-09,-2.0864e-13,7133.44,-2.22829], Tmin=(776.878,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.6206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cdd-CdsCds) + radical(CCOJ)"""),
)

species(
    label = '[O]C=[C]C[C]([O])O[O](24693)',
    structure = SMILES('[O]C=[C]C[C]([O])O[O]'),
    E0 = (448.393,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,360,370,350,3010,987.5,1337.5,450,1655,2750,2850,1437.5,1250,1305,750,350,336.58,338.005,338.377,2959.31],'cm^-1')),
        HinderedRotor(inertia=(0.136326,'amu*angstrom^2'), symmetry=1, barrier=(11.1718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137056,'amu*angstrom^2'), symmetry=1, barrier=(11.1601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.397597,'amu*angstrom^2'), symmetry=1, barrier=(32.6572,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696516,0.0782745,-0.000121345,1.00163e-07,-3.28838e-11,54042.8,32.612], Tmin=(100,'K'), Tmax=(790.928,'K')), NASAPolynomial(coeffs=[10.6896,0.0241801,-1.20113e-05,2.32255e-09,-1.61416e-13,52573.2,-12.552], Tmin=(790.928,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(CCOJ) + radical(C=COJ) + radical(Cds_S) + radical(Cs_P)"""),
)

species(
    label = '[O]C=C[C]C([O])O[O](24694)',
    structure = SMILES('[O]C=C[C]C([O])O[O]'),
    E0 = (453.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,370.169,370.169,370.169,370.169,370.169],'cm^-1')),
        HinderedRotor(inertia=(0.129942,'amu*angstrom^2'), symmetry=1, barrier=(12.635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129942,'amu*angstrom^2'), symmetry=1, barrier=(12.635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.417675,'amu*angstrom^2'), symmetry=1, barrier=(40.613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81159,0.0738024,-9.70609e-05,6.57203e-08,-1.76963e-11,54657.5,30.9537], Tmin=(100,'K'), Tmax=(907.057,'K')), NASAPolynomial(coeffs=[12.6924,0.0214091,-1.04178e-05,2.03915e-09,-1.44596e-13,52502.2,-25.2055], Tmin=(907.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(CCJ2_triplet) + radical(ROOJ)"""),
)

species(
    label = '[O][C]=C[CH]C([O])O[O](24695)',
    structure = SMILES('[O][C]=C[CH]C([O])O[O]'),
    E0 = (361.965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,3010,987.5,1337.5,450,1655,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,316.453,316.683,316.843,1765.28],'cm^-1')),
        HinderedRotor(inertia=(0.193722,'amu*angstrom^2'), symmetry=1, barrier=(13.5683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.639441,'amu*angstrom^2'), symmetry=1, barrier=(44.8053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636763,'amu*angstrom^2'), symmetry=1, barrier=(44.8096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00544,0.0739729,-9.8914e-05,6.15472e-08,-9.46095e-12,43634.5,30.7192], Tmin=(100,'K'), Tmax=(591.069,'K')), NASAPolynomial(coeffs=[9.74831,0.0276848,-1.41277e-05,2.77933e-09,-1.95869e-13,42376,-8.76591], Tmin=(591.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CCJCO) + radical(C=CJO) + radical(CCOJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C[CH][C]([O])O[O](24696)',
    structure = SMILES('[O]C=C[CH][C]([O])O[O]'),
    E0 = (327.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,400.523,400.533,400.547,400.547],'cm^-1')),
        HinderedRotor(inertia=(0.167871,'amu*angstrom^2'), symmetry=1, barrier=(19.1125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289038,'amu*angstrom^2'), symmetry=1, barrier=(32.9126,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.552012,'amu*angstrom^2'), symmetry=1, barrier=(62.8576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.733331,0.0748465,-9.37015e-05,5.94787e-08,-1.49956e-11,39500.4,29.2158], Tmin=(100,'K'), Tmax=(966.631,'K')), NASAPolynomial(coeffs=[13.6968,0.0212034,-1.04602e-05,2.06956e-09,-1.4805e-13,36994.2,-32.8855], Tmin=(966.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_P) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C[C][CH]C(=O)O[O](24697)',
    structure = SMILES('[O]C[C][CH]C(=O)O[O]'),
    E0 = (413.433,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21145,0.0661372,-8.58173e-05,6.33607e-08,-1.9331e-11,49820.5,31.9202], Tmin=(100,'K'), Tmax=(792.802,'K')), NASAPolynomial(coeffs=[8.77363,0.0279819,-1.36249e-05,2.65252e-09,-1.86924e-13,48621.5,-2.807], Tmin=(792.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.433,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + radical(CCJ2_triplet) + radical(CCOJ) + radical(CCJCO) + radical(C(=O)OOJ)"""),
)

species(
    label = 'O2(S)(5486)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121571,5.31618e-06,-4.89443e-09,1.45845e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.56,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69971e-07,1.51275e-10,-1.08782e-14,10302.3,6.16754], Tmin=(1074.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
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
    label = '[O]C=C=CC1OO1(24698)',
    structure = SMILES('[O]C=C=CC1OO1'),
    E0 = (113.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37795,0.0411526,8.35655e-06,-5.82479e-08,3.07861e-11,13763.8,22.1193], Tmin=(100,'K'), Tmax=(905.925,'K')), NASAPolynomial(coeffs=[20.7725,-0.000223836,3.58555e-06,-8.09927e-10,5.3729e-14,8433.66,-79.5557], Tmin=(905.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.525,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(dioxirane) + radical(C=COJ)"""),
)

species(
    label = '[O]C1C=C=COO1(24699)',
    structure = SMILES('[O]C1C=C=COO1'),
    E0 = (177.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.380677,0.10944,-0.00019765,1.82181e-07,-6.44316e-11,21491.9,10.8795], Tmin=(100,'K'), Tmax=(836.117,'K')), NASAPolynomial(coeffs=[9.79788,0.0354122,-1.93948e-05,3.81425e-09,-2.64591e-13,20675.4,-31.1089], Tmin=(836.117,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(six-inringtwodouble-12) + radical(CCOJ)"""),
)

species(
    label = '[O]O(16)',
    structure = SMILES('[O]O'),
    E0 = (-8.19602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1036.72,2034.11,2034.11],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.04595,-0.00173474,1.0377e-05,-1.02207e-08,3.3493e-12,-986.755,4.63579], Tmin=(100,'K'), Tmax=(932.129,'K')), NASAPolynomial(coeffs=[3.21022,0.00367946,-1.27704e-06,2.18051e-10,-1.46343e-14,-910.359,8.18305], Tmin=(932.129,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.19602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsH) + group(O2s-OsH) + radical(HOOJ)"""),
)

species(
    label = '[O]C=C=C=C[O](24700)',
    structure = SMILES('[O]C=C=C=C[O]'),
    E0 = (188.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49248,0.0406836,-6.88737e-06,-3.893e-08,2.32066e-11,22792.8,17.926], Tmin=(100,'K'), Tmax=(911.655,'K')), NASAPolynomial(coeffs=[20.6682,-0.00488084,4.61894e-06,-9.35092e-10,6.07591e-14,17693.6,-81.6029], Tmin=(911.655,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CC1OOO1(24701)',
    structure = SMILES('[O]C=C=CC1OOO1'),
    E0 = (154.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28536,0.0371853,3.66733e-05,-8.71592e-08,3.88522e-11,18728.8,25.5458], Tmin=(100,'K'), Tmax=(956.366,'K')), NASAPolynomial(coeffs=[22.4165,0.00491193,-7.09076e-07,2.43443e-10,-3.1146e-14,12121.1,-88.8719], Tmin=(956.366,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclobutane) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C=CC([O])[O](24702)',
    structure = SMILES('[O]C=C=CC([O])[O]'),
    E0 = (170.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,1517.23],'cm^-1')),
        HinderedRotor(inertia=(0.980577,'amu*angstrom^2'), symmetry=1, barrier=(22.5454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11534,0.0503278,-3.6232e-05,-3.81381e-08,6.02561e-11,20592,23.0944], Tmin=(100,'K'), Tmax=(468.994,'K')), NASAPolynomial(coeffs=[6.67055,0.0274256,-1.3992e-05,2.75677e-09,-1.94306e-13,19989.4,2.69725], Tmin=(468.994,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CCOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]C=C=C[CH]O[O](24703)',
    structure = SMILES('[O]C=C=C[CH]O[O]'),
    E0 = (267.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.45998,'amu*angstrom^2'), symmetry=1, barrier=(33.5678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.45833,'amu*angstrom^2'), symmetry=1, barrier=(33.5299,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20909,0.0504798,-2.83308e-05,-1.0559e-08,1.06472e-11,32228.1,25.1849], Tmin=(100,'K'), Tmax=(948.007,'K')), NASAPolynomial(coeffs=[17.8289,0.00697677,-1.62096e-06,2.80405e-10,-2.30939e-14,27880.7,-60.4184], Tmin=(948.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CCJO) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C=CC([O])O[O](19787)',
    structure = SMILES('[CH]=C=CC([O])O[O]'),
    E0 = (390.335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,492.5,1135,1000,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180],'cm^-1')),
        HinderedRotor(inertia=(0.442811,'amu*angstrom^2'), symmetry=1, barrier=(10.1811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.438226,'amu*angstrom^2'), symmetry=1, barrier=(10.0757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28906,0.0695021,-0.000123677,1.16366e-07,-4.13208e-11,47034.7,26.6256], Tmin=(100,'K'), Tmax=(871.807,'K')), NASAPolynomial(coeffs=[5.38325,0.0281482,-1.36936e-05,2.56749e-09,-1.72855e-13,47178.5,12.3541], Tmin=(871.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCOJ) + radical(C=C=CJ)"""),
)

species(
    label = '[O]OC([O])C#CC=O(24704)',
    structure = SMILES('[O]OC([O])C#CC=O'),
    E0 = (131.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2100,2250,500,550,2782.5,750,1395,475,1775,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.209408,'amu*angstrom^2'), symmetry=1, barrier=(4.8147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.209262,'amu*angstrom^2'), symmetry=1, barrier=(4.81134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.42733,'amu*angstrom^2'), symmetry=1, barrier=(55.8091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.056,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13967,0.0743234,-0.000136622,1.31554e-07,-4.76652e-11,15880.1,28.1186], Tmin=(100,'K'), Tmax=(861.407,'K')), NASAPolynomial(coeffs=[4.9093,0.0308592,-1.57315e-05,3.0083e-09,-2.04961e-13,16193.8,16.0852], Tmin=(861.407,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[O]OC([O])[C]=CC=O(24705)',
    structure = SMILES('[O]OC([O])[C]=CC=O'),
    E0 = (209.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,492.5,1135,1000,2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.260174,'amu*angstrom^2'), symmetry=1, barrier=(5.98191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.258844,'amu*angstrom^2'), symmetry=1, barrier=(5.95133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261337,'amu*angstrom^2'), symmetry=1, barrier=(6.00865,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889956,0.0846224,-0.00016869,1.74528e-07,-6.7165e-11,25309.4,30.1164], Tmin=(100,'K'), Tmax=(837.247,'K')), NASAPolynomial(coeffs=[2.09698,0.0410662,-2.29528e-05,4.5734e-09,-3.19683e-13,26431.8,32.4174], Tmin=(837.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(Cds_S) + radical(ROOJ)"""),
)

species(
    label = '[O]OC([O])[CH]C=C=O(24706)',
    structure = SMILES('[O]OC([O])[CH]C=C=O'),
    E0 = (104.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.739082,0.0808729,-0.000136471,1.25366e-07,-4.49985e-11,12663.7,27.9591], Tmin=(100,'K'), Tmax=(819.906,'K')), NASAPolynomial(coeffs=[7.33167,0.0320749,-1.67617e-05,3.28438e-09,-2.28888e-13,12141.8,0.872549], Tmin=(819.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCOJ) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C=C[CH]C(=O)O[O](24707)',
    structure = SMILES('[O]C=C[CH]C(=O)O[O]'),
    E0 = (-117.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82922,0.0487252,1.0249e-05,-6.14174e-08,2.95744e-11,-13940.1,25.8124], Tmin=(100,'K'), Tmax=(974.331,'K')), NASAPolynomial(coeffs=[23.5887,0.0060546,-2.2138e-06,5.8599e-10,-5.60264e-14,-20784.8,-95.7628], Tmin=(974.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C(=O)OOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OC1C=[C][CH]OO1(24626)',
    structure = SMILES('[O]OC1C=[C][CH]OO1'),
    E0 = (340.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47847,0.0406039,6.48473e-06,-4.29192e-08,2.02003e-11,41018.6,22.8557], Tmin=(100,'K'), Tmax=(987.386,'K')), NASAPolynomial(coeffs=[16.9811,0.0130272,-5.13616e-06,1.059e-09,-8.30873e-14,36240,-60.4343], Tmin=(987.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.183,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(36dihydro12dioxin) + radical(C=CCJO) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C1C=[C][CH]OOO1(24708)',
    structure = SMILES('[O]C1C=[C][CH]OOO1'),
    E0 = (395.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99574,-0.00926331,0.000216504,-3.14853e-07,1.31696e-10,47648.7,28.0509], Tmin=(100,'K'), Tmax=(910.121,'K')), NASAPolynomial(coeffs=[38.2792,-0.0297809,2.13139e-05,-4.12701e-09,2.65084e-13,35289.5,-175.194], Tmin=(910.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(CCOJ) + radical(Cds_S) + radical(C=CCJO)"""),
)

species(
    label = '[O]OC(=O)C=CC=O(24709)',
    structure = SMILES('[O]OC(=O)C=CC=O'),
    E0 = (-185.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45187,0.0624379,-8.02225e-05,5.90051e-08,-1.8352e-11,-22177,27.2051], Tmin=(100,'K'), Tmax=(769.446,'K')), NASAPolynomial(coeffs=[7.80136,0.0294291,-1.58721e-05,3.24914e-09,-2.36033e-13,-23154.1,-1.76327], Tmin=(769.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-185.105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cd-Cd(CO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)O2s) + group(Cds-O2d(Cds-Cds)H) + radical(C(=O)OOJ)"""),
)

species(
    label = '[O]C1C=C(C=O)O1(24710)',
    structure = SMILES('[O]C1C=C(C=O)O1'),
    E0 = (-97.6506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47786,0.0465306,-3.08055e-05,4.51029e-09,1.47531e-12,-11646.4,19.2707], Tmin=(100,'K'), Tmax=(1171.18,'K')), NASAPolynomial(coeffs=[15.2842,0.0137127,-7.13408e-06,1.48708e-09,-1.10256e-13,-15863.5,-53.7157], Tmin=(1171.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.6506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(CCOJ)"""),
)

species(
    label = '[O]OC1C=C(C=O)O1(24577)',
    structure = SMILES('[O]OC1C=C(C=O)O1'),
    E0 = (-99.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.670246,0.0630976,-5.65125e-05,2.06229e-08,-2.03165e-12,-11880.2,22.8243], Tmin=(100,'K'), Tmax=(1126.65,'K')), NASAPolynomial(coeffs=[19.1214,0.011184,-5.49492e-06,1.14427e-09,-8.58294e-14,-16900.6,-72.2209], Tmin=(1126.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(ROOJ)"""),
)

species(
    label = '[O]C1C=C(C=O)OO1(24711)',
    structure = SMILES('[O]C1C=C(C=O)OO1'),
    E0 = (-96.819,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.064,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51656,0.0449091,-3.03125e-05,8.99555e-09,-1.09381e-12,-11603.4,20.9472], Tmin=(100,'K'), Tmax=(1680.04,'K')), NASAPolynomial(coeffs=[9.44248,0.0284192,-1.55896e-05,3.15328e-09,-2.24448e-13,-13930.5,-16.0597], Tmin=(1680.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-96.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + ring(12dioxolene) + radical(CCOJ)"""),
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
    label = '[C]=CC([O])O[O](8561)',
    structure = SMILES('[C]=CC([O])O[O]'),
    E0 = (653.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,180,2052.46],'cm^-1')),
        HinderedRotor(inertia=(0.376272,'amu*angstrom^2'), symmetry=1, barrier=(8.65124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.377304,'amu*angstrom^2'), symmetry=1, barrier=(8.67497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86484,0.0594853,-0.000122945,1.27074e-07,-4.79905e-11,78648.6,23.4488], Tmin=(100,'K'), Tmax=(861.942,'K')), NASAPolynomial(coeffs=[2.24539,0.0275562,-1.48883e-05,2.89866e-09,-1.98788e-13,79703.4,28.1691], Tmin=(861.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CdCdJ2_triplet) + radical(ROOJ)"""),
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
    E0 = (168.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (252.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (290.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (255.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (427.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (168.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (168.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (293.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (444.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (333.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (322.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (314.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (561.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (259.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (259.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (241.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (168.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (498.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (332.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (618.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (518.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (355.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (294.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (294.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (251.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (238.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (280.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (231.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (193.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (193.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (193.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (474.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (475.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (384.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (390.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (421.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (168.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (358.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (431.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (347.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (176.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (418.925,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (673.923,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (797.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (352.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (329.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (406.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (331.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (345.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (340.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (395.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (257.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (251.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (176.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (175.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (721.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC=O(5472)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=[C]C1OC1O[O](24581)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(83.6851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=[C]C1OOC1[O](24670)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', '[O]C=C=CC(=O)O[O](24671)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(25.1243,'m^3/(mol*s)'), n=1.86, Ea=(32.426,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO-DeNd_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[O]OC([O])C=C=C=O(24672)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]OC=O(5472)', '[CH]=C=C[O](8556)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(94.1495,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-NdH_O;CdsJ=Cdd]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond
Ea raised from 92.1 to 94.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['O2(2)', '[O]C=C=CC=O(22476)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.5e+07,'cm^3/(mol*s)'), n=2.16, Ea=(185.246,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CdH_O;YJ] for rate rule [CO-CdH_O;O2b]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 181.8 to 185.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=[C]C=C(O)O[O](24673)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.66008e+10,'s^-1'), n=0.776642, Ea=(125.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;Cs_H_out_noH] + [R2H_S;O_rad_out;Cs_H_out] for rate rule [R2H_S;O_rad_out;Cs_H_out_noH]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC([O])C=C=[C]O(24674)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=[C][CH]C(=O)OO(24675)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.83109e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C=C=[C]C(O)O[O](24676)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=C=[C]C([O])OO(24677)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 337 used for R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC([O])C#C[CH]O(24678)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC(O)[CH][C]=C=O(24679)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC(=O)[CH][C]=CO(24680)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(246072,'s^-1'), n=1.73305, Ea=(91.0472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;O_rad_out;XH_out] for rate rule [R5H_SMMS;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C([CH][C]=C=O)OO(24681)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.33753e+06,'s^-1'), n=1.02312, Ea=(72.6006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6H;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[O]C=[C]C=C[O](23191)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(18.8989,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 18.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][CH]O[O](8201)', '[CH]=C=C[O](8556)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[O]C=[C][CH]C(=O)O[O](24682)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[O]C=C=[C]C([O])O[O](24683)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(8)', '[O]OC([O])[CH][C]=C=O(24684)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC([O])C=C1[CH]O1(24685)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=C1[CH]C(O[O])O1(24686)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC([O])C1[C]=CO1(24687)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=C1[CH]C([O])OO1(24688)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.66666e+10,'s^-1'), n=0.302034, Ea=(82.5645,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC1C=[C]C([O])O1(24641)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.01304e+12,'s^-1'), n=-0.3725, Ea=(69.6427,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C1[C]=CC([O])OO1(24689)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.26053e+10,'s^-1'), n=0.2505, Ea=(112.349,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R6_linear;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=C=CC(=O)OO(22529)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC(=O)C=C=CO(24690)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC(O)C=C=C=O(24691)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C(C=C=C=O)OO(24692)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C=[C]C[C]([O])O[O](24693)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.64935e+10,'s^-1'), n=0.2551, Ea=(25.8154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C=C[C]C([O])O[O](24694)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O][C]=C[CH]C([O])O[O](24695)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C=C[CH][C]([O])O[O](24696)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C[C][CH]C(=O)O[O](24697)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['O2(S)(5486)', '[O]C=C=CC=O(22476)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['O(T)(63)', '[O]C=C=CC1OO1(24698)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['O(T)(63)', '[O]C1C=C=COO1(24699)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(9.23936e+08,'s^-1'), n=0.739184, Ea=(263.098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R5OO;Y_rad_intra;OO] for rate rule [R5OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]O(16)', '[O]C=C=C=C[O](24700)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=C=CC1OOO1(24701)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['O(T)(63)', '[O]C=C=CC([O])[O](24702)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['O(T)(63)', '[O]C=C=C[CH]O[O](24703)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['O(T)(63)', '[CH]=C=CC([O])O[O](19787)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['H(8)', '[O]OC([O])C#CC=O(24704)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(3249.46,'m^3/(mol*s)'), n=1.38433, Ea=(9.80868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-De;HJ] for rate rule [Ct-Cs_Ct-CO;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[O][CH]O[O](8201)', 'C#CC=O(21959)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.106247,'m^3/(mol*s)'), n=2.32278, Ea=(16.475,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-De;CJ] for rate rule [Ct-H_Ct-CO;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[O]OC([O])[C]=CC=O(24705)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC([O])[CH]C=C=O(24706)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C=C[CH]C(=O)O[O](24707)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(8.2826e+06,'s^-1'), n=1.67955, Ea=(176.633,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleDe;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC1C=[C][CH]OO1(24626)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(171.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic
Ea raised from 167.0 to 171.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction51',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C1C=[C][CH]OOO1(24708)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(226.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 216.4 to 226.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC(=O)C=CC=O(24709)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['O(T)(63)', '[O]C1C=C(C=O)O1(24710)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO] for rate rule [R3OO_DS;Cd_rad_in/OneDe;OOJ]
Euclidian distance = 3.31662479036
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]OC1C=C(C=O)O1(24577)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[O]C=C=CC([O])O[O](22518)'],
    products = ['[O]C1C=C(C=O)OO1(24711)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;O_rad;CdsinglepriDe_rad_out]
Euclidian distance = 2.44948974278
family: Birad_recombination"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]=O(373)', '[C]=CC([O])O[O](8561)'],
    products = ['[O]C=C=CC([O])O[O](22518)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4819',
    isomers = [
        '[O]C=C=CC([O])O[O](22518)',
    ],
    reactants = [
        ('[O]OC=O(5472)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4819',
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

