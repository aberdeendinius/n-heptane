species(
    label = '[CH]=C(C=O)C=C=C[O](22652)',
    structure = SMILES('[CH]=C(C=O)C=C=C[O]'),
    E0 = (300.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180],'cm^-1')),
        HinderedRotor(inertia=(1.26441,'amu*angstrom^2'), symmetry=1, barrier=(29.0713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26476,'amu*angstrom^2'), symmetry=1, barrier=(29.0794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0723359,0.0739509,-6.32921e-05,1.35797e-08,3.91554e-12,36271.8,25.8627], Tmin=(100,'K'), Tmax=(977.162,'K')), NASAPolynomial(coeffs=[23.1468,0.00679682,-2.11532e-06,4.34057e-10,-3.63074e-14,30458.9,-91.5945], Tmin=(977.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(Cds_P) + radical(C=COJ)"""),
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
    label = '[O]C=C=CC1=CC1[O](25078)',
    structure = SMILES('[O]C=C=CC1=CC1[O]'),
    E0 = (450.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.38113,0.0626772,-2.57468e-05,-3.01562e-08,2.10877e-11,54344.7,25.4951], Tmin=(100,'K'), Tmax=(935.574,'K')), NASAPolynomial(coeffs=[24.2951,0.0032582,8.59766e-07,-1.90366e-10,6.81124e-15,47995.8,-98.2999], Tmin=(935.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclopropene) + radical(C=COJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C1C=C=COC1[O](25079)',
    structure = SMILES('[CH]=C1C=C=COC1[O]'),
    E0 = (339.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10341,0.0413977,3.46596e-05,-8.31361e-08,3.61557e-11,40942.5,17.9488], Tmin=(100,'K'), Tmax=(975.83,'K')), NASAPolynomial(coeffs=[21.4792,0.0124345,-4.68449e-06,1.03789e-09,-8.75239e-14,34368.2,-93.1641], Tmin=(975.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(CCOJ) + radical(Cds_P)"""),
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
    label = '[CH]=C(C=O)C=C=C=O(25080)',
    structure = SMILES('[CH]=C(C=O)C=C=C=O'),
    E0 = (344.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,540,610,2055,350,440,435,1725,2120,512.5,787.5,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.000808762,'amu*angstrom^2'), symmetry=1, barrier=(9.18273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000807907,'amu*angstrom^2'), symmetry=1, barrier=(9.17302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33476,0.0505303,-4.51588e-05,1.53329e-08,-8.13851e-13,41525.8,8.5368], Tmin=(100,'K'), Tmax=(1071.73,'K')), NASAPolynomial(coeffs=[16.3231,0.00764125,-3.39799e-06,7.18652e-10,-5.54303e-14,37563.5,-68.3083], Tmin=(1071.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P)"""),
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
    label = 'C#CC=C=C[O](23376)',
    structure = SMILES('C#CC=C=C[O]'),
    E0 = (343.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.6143,'amu*angstrom^2'), symmetry=1, barrier=(37.1158,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58392,0.0406745,-9.28843e-06,-2.95441e-08,1.76387e-11,41466.8,17.6904], Tmin=(100,'K'), Tmax=(933.363,'K')), NASAPolynomial(coeffs=[17.8115,0.00285504,5.05668e-07,-1.22779e-10,4.01942e-15,37055.7,-66.882], Tmin=(933.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsCtH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(C=O)=CC#CO(25081)',
    structure = SMILES('[CH]C(C=O)=CC#CO'),
    E0 = (333.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3615,1277.5,1000,2100,2250,500,550,350,440,435,1725,187.761,187.781,187.781,187.818],'cm^-1')),
        HinderedRotor(inertia=(1.90886,'amu*angstrom^2'), symmetry=1, barrier=(47.7565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9089,'amu*angstrom^2'), symmetry=1, barrier=(47.7569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90881,'amu*angstrom^2'), symmetry=1, barrier=(47.7567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90834,'amu*angstrom^2'), symmetry=1, barrier=(47.7566,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40503,0.0618128,-4.85083e-05,1.89043e-08,-3.10193e-12,40151,25.0744], Tmin=(100,'K'), Tmax=(1359.64,'K')), NASAPolynomial(coeffs=[11.1976,0.0330034,-1.67248e-05,3.32011e-09,-2.36426e-13,37488.1,-25.1778], Tmin=(1359.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CtH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtOs) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(C#CC=O)=C[O](25082)',
    structure = SMILES('[CH2]C(C#CC=O)=C[O]'),
    E0 = (184.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769888,0.0601164,-3.90921e-05,-1.0217e-09,6.87586e-12,22311.5,26.4336], Tmin=(100,'K'), Tmax=(979.556,'K')), NASAPolynomial(coeffs=[18.1557,0.0133937,-4.7134e-06,8.76676e-10,-6.45865e-14,17740.9,-63.028], Tmin=(979.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([C]=O)C=C=C[O](25083)',
    structure = SMILES('C=C([C]=O)C=C=C[O]'),
    E0 = (216.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0940752,0.0735734,-6.29616e-05,1.3604e-08,3.836e-12,26229.2,25.943], Tmin=(100,'K'), Tmax=(977.479,'K')), NASAPolynomial(coeffs=[22.9932,0.00691018,-2.16293e-06,4.41948e-10,-3.67714e-14,20460.5,-90.6192], Tmin=(977.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=C(C)CJ=O) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=C=C=CO)C=O(25084)',
    structure = SMILES('[CH]C(=C=C=CO)C=O'),
    E0 = (316.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,540,563.333,586.667,610,1970,2140,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.87176,'amu*angstrom^2'), symmetry=1, barrier=(43.0355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86936,'amu*angstrom^2'), symmetry=1, barrier=(42.9802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.88392,'amu*angstrom^2'), symmetry=1, barrier=(43.3151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.285541,0.0843414,-8.43796e-05,3.68176e-08,-4.9799e-12,38201.2,26.6316], Tmin=(100,'K'), Tmax=(1019.65,'K')), NASAPolynomial(coeffs=[21.7668,0.0142554,-5.43733e-06,1.00051e-09,-7.12958e-14,32850.3,-84.374], Tmin=(1019.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(C=C=C=O)=C[O](25085)',
    structure = SMILES('[CH2]C(C=C=C=O)=C[O]'),
    E0 = (257.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26663,0.043574,2.32003e-06,-5.06589e-08,2.71977e-11,31125.3,9.28835], Tmin=(100,'K'), Tmax=(924.146,'K')), NASAPolynomial(coeffs=[21.1505,0.000467653,2.56186e-06,-5.34626e-10,3.13027e-14,25615.8,-94.9961], Tmin=(924.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]C(=C=O)C=C=CO(25086)',
    structure = SMILES('[CH]C(=C=O)C=C=CO'),
    E0 = (302.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.140293,0.0812805,-9.1823e-05,5.29533e-08,-1.19725e-11,36529.4,27.0458], Tmin=(100,'K'), Tmax=(1085.61,'K')), NASAPolynomial(coeffs=[16.6441,0.0204711,-7.80153e-06,1.356e-09,-9.03742e-14,32946,-53.931], Tmin=(1085.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + group(Cd-Cd(CCO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C=C[C]=C[O](23385)',
    structure = SMILES('[CH]=C=C[C]=C[O]'),
    E0 = (520.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.70898,'amu*angstrom^2'), symmetry=1, barrier=(39.2927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.518732,0.0587054,-6.62874e-05,3.49603e-08,-6.66837e-12,62759.5,22.8411], Tmin=(100,'K'), Tmax=(1556.85,'K')), NASAPolynomial(coeffs=[16.6546,0.0022419,2.57223e-06,-7.17553e-10,5.49875e-14,59553.8,-56.3074], Tmin=(1556.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(520.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJC=C) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C(C#CC=O)=C[O](25087)',
    structure = SMILES('[CH]C(C#CC=O)=C[O]'),
    E0 = (403.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,350,440,435,1725,2100,2250,500,550,3010,987.5,1337.5,450,1655,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10927,'amu*angstrom^2'), symmetry=1, barrier=(48.4962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11218,'amu*angstrom^2'), symmetry=1, barrier=(48.5631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.11033,'amu*angstrom^2'), symmetry=1, barrier=(48.5207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744596,0.0628215,-4.81785e-05,1.28073e-08,7.54744e-13,48672.1,27.3533], Tmin=(100,'K'), Tmax=(1027.8,'K')), NASAPolynomial(coeffs=[16.0015,0.019213,-7.54806e-06,1.38001e-09,-9.65529e-14,44703.1,-50.7226], Tmin=(1027.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(403.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCtCs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C=O)C=C=C[O](25088)',
    structure = SMILES('[CH]C(=C=O)C=C=C[O]'),
    E0 = (444.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,350,440,435,1725,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09904,'amu*angstrom^2'), symmetry=1, barrier=(48.261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10081,'amu*angstrom^2'), symmetry=1, barrier=(48.3017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.780547,0.0715302,-7.83266e-05,4.55802e-08,-1.07075e-11,53516.6,26.164], Tmin=(100,'K'), Tmax=(1029.35,'K')), NASAPolynomial(coeffs=[12.6406,0.0254442,-1.11708e-05,2.08764e-09,-1.44686e-13,51075,-31.3972], Tmin=(1029.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + group(Cd-Cd(CCO)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C=O)=C[C]=C=O(25089)',
    structure = SMILES('[CH]C(C=O)=C[C]=C=O'),
    E0 = (470.109,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2120,512.5,787.5,350,440,435,1725,3010,987.5,1337.5,450,1655,258.776,261.743,261.92,262.682],'cm^-1')),
        HinderedRotor(inertia=(1.02323,'amu*angstrom^2'), symmetry=1, barrier=(50.7226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05252,'amu*angstrom^2'), symmetry=1, barrier=(50.741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03888,'amu*angstrom^2'), symmetry=1, barrier=(50.7163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.894163,0.0782705,-0.000126914,1.21409e-07,-4.58926e-11,56643.2,27.4664], Tmin=(100,'K'), Tmax=(807.97,'K')), NASAPolynomial(coeffs=[4.41233,0.0411568,-2.14457e-05,4.21408e-09,-2.95064e-13,56717.6,15.222], Tmin=(807.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C=O)C=C1[CH]O1(25090)',
    structure = SMILES('[CH]=C(C=O)C=C1[CH]O1'),
    E0 = (326.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.478975,0.0529282,1.3345e-05,-7.81104e-08,3.95286e-11,39417.6,22.9844], Tmin=(100,'K'), Tmax=(937.804,'K')), NASAPolynomial(coeffs=[28.2669,-0.00339077,3.9312e-06,-6.89226e-10,3.47027e-14,31470.3,-123.876], Tmin=(937.804,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_P) + radical(C=CCJO)"""),
)

species(
    label = 'O=C[C]1[CH]C(C=O)=C1(25091)',
    structure = SMILES('O=C[C]1[CH]C(C=O)=C1'),
    E0 = (219.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9209,0.03241,2.1218e-05,-4.68262e-08,1.83405e-11,26520.9,23.4331], Tmin=(100,'K'), Tmax=(1050.7,'K')), NASAPolynomial(coeffs=[12.7301,0.0230084,-1.06859e-05,2.17612e-09,-1.61877e-13,22496.9,-37.5895], Tmin=(1050.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.789,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(CCJCC=O) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[O]C=C=CC1[CH]OC=1(25092)',
    structure = SMILES('[O]C=C=CC1[CH]OC=1'),
    E0 = (322.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0330309,0.0590638,1.52603e-05,-9.67405e-08,5.13332e-11,38921.3,22.6259], Tmin=(100,'K'), Tmax=(907.925,'K')), NASAPolynomial(coeffs=[33.6091,-0.0135669,1.08606e-05,-2.17027e-09,1.42268e-13,29721,-153.206], Tmin=(907.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclobutene) + radical(C=COJ) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH]=C(C=O)C1[C]=CO1(25093)',
    structure = SMILES('[CH]=C(C=O)C1[C]=CO1'),
    E0 = (416.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.709471,0.0556758,-1.59656e-05,-3.03518e-08,1.78031e-11,50272.5,25.6462], Tmin=(100,'K'), Tmax=(988.709,'K')), NASAPolynomial(coeffs=[21.8259,0.0088757,-3.57075e-06,8.08234e-10,-6.80852e-14,44208.8,-85.5371], Tmin=(988.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.88,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_P) + radical(Cds_S)"""),
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
    label = '[CH]=C1[CH]OOC=C=C1(25095)',
    structure = SMILES('[CH]=C1[CH]OOC=C=C1'),
    E0 = (624.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16928,0.0424995,2.41498e-05,-7.08933e-08,3.21368e-11,75249.3,19.2385], Tmin=(100,'K'), Tmax=(963.062,'K')), NASAPolynomial(coeffs=[20.2661,0.0117744,-3.67841e-06,7.61333e-10,-6.40543e-14,69317.6,-83.8732], Tmin=(963.062,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(1_2_cycloheptadiene) + radical(C=CCJO) + radical(Cds_P)"""),
)

species(
    label = 'C=C(C=O)C=C=C=O(25096)',
    structure = SMILES('C=C(C=O)C=C=C=O'),
    E0 = (97.3129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42367,0.0468407,-2.79792e-05,-2.62546e-09,5.29891e-12,11805.4,7.81386], Tmin=(100,'K'), Tmax=(1034.71,'K')), NASAPolynomial(coeffs=[15.9364,0.0107139,-4.56659e-06,9.48577e-10,-7.28431e-14,7732.69,-67.8643], Tmin=(1034.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.3129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]C(=[C]C=C[O])C=O(25097)',
    structure = SMILES('[CH]C(=[C]C=C[O])C=O'),
    E0 = (429.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0379,'amu*angstrom^2'), symmetry=1, barrier=(46.8554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03623,'amu*angstrom^2'), symmetry=1, barrier=(46.817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03697,'amu*angstrom^2'), symmetry=1, barrier=(46.834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.212382,0.0766644,-7.52635e-05,3.70888e-08,-7.21197e-12,51791.6,27.9983], Tmin=(100,'K'), Tmax=(1249.34,'K')), NASAPolynomial(coeffs=[17.5151,0.0212659,-8.74933e-06,1.59539e-09,-1.09463e-13,47468.3,-59.3288], Tmin=(1249.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C=C[C]=O)=C[O](25098)',
    structure = SMILES('[CH]C(C=C[C]=O)=C[O]'),
    E0 = (396.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.92216,'amu*angstrom^2'), symmetry=1, barrier=(44.1943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9312,'amu*angstrom^2'), symmetry=1, barrier=(44.4021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.9272,'amu*angstrom^2'), symmetry=1, barrier=(44.3101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.185458,0.0848095,-9.11462e-05,4.76902e-08,-9.68896e-12,47883.1,28.6008], Tmin=(100,'K'), Tmax=(1208.01,'K')), NASAPolynomial(coeffs=[20.2503,0.0171402,-7.11833e-06,1.31637e-09,-9.15774e-14,42945.9,-73.8514], Tmin=(1208.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C=C=C[O])C[O](25099)',
    structure = SMILES('[CH]C(=C=C=C[O])C[O]'),
    E0 = (624.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,563.333,586.667,610,1970,2140,350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.1722,'amu*angstrom^2'), symmetry=1, barrier=(49.9431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.17284,'amu*angstrom^2'), symmetry=1, barrier=(49.9579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.745547,0.0754257,-8.54836e-05,5.47417e-08,-1.45399e-11,75181.8,28.3306], Tmin=(100,'K'), Tmax=(904.267,'K')), NASAPolynomial(coeffs=[10.3325,0.0330173,-1.51352e-05,2.87677e-09,-2.00694e-13,73448,-16.9561], Tmin=(904.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(CCOJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C=O)C[C]=C[O](25100)',
    structure = SMILES('[CH]C(=C=O)C[C]=C[O]'),
    E0 = (513.026,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,350,440,435,1725,3010,987.5,1337.5,450,1655,252.307,252.36,252.362,252.366,252.366,252.42],'cm^-1')),
        HinderedRotor(inertia=(1.11733,'amu*angstrom^2'), symmetry=1, barrier=(50.5153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11785,'amu*angstrom^2'), symmetry=1, barrier=(50.5123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11687,'amu*angstrom^2'), symmetry=1, barrier=(50.5111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.531941,0.0785339,-9.26157e-05,5.91819e-08,-1.52777e-11,61825.8,28.4492], Tmin=(100,'K'), Tmax=(940.618,'K')), NASAPolynomial(coeffs=[12.3683,0.0282001,-1.23495e-05,2.29349e-09,-1.57884e-13,59599.1,-27.9299], Tmin=(940.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(513.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + other(ketene_1C-C_1C-H) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=C=C=C[O])[CH]O(25101)',
    structure = SMILES('[CH]C(=C=C=C[O])[CH]O'),
    E0 = (515.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,540,563.333,586.667,610,1970,2140,3010,987.5,1337.5,450,1655,3615,1277.5,1000,350,440,435,1725,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06928,'amu*angstrom^2'), symmetry=1, barrier=(47.5768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07063,'amu*angstrom^2'), symmetry=1, barrier=(47.6078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06805,'amu*angstrom^2'), symmetry=1, barrier=(47.5486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199283,0.0719419,-5.10236e-05,3.49788e-09,6.87954e-12,62176.6,30.0136], Tmin=(100,'K'), Tmax=(959.758,'K')), NASAPolynomial(coeffs=[20.1643,0.0156749,-5.19075e-06,9.09571e-10,-6.48574e-14,57103.4,-71.9503], Tmin=(959.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(AllylJ2_triplet) + radical(C=CCJO) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=C=O)[CH]C=C[O](25102)',
    structure = SMILES('[CH]C(=C=O)[CH]C=C[O]'),
    E0 = (376.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2120,512.5,787.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,257.551,257.552,257.553,257.553,257.553,257.554],'cm^-1')),
        HinderedRotor(inertia=(1.05004,'amu*angstrom^2'), symmetry=1, barrier=(49.4274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05005,'amu*angstrom^2'), symmetry=1, barrier=(49.4274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05005,'amu*angstrom^2'), symmetry=1, barrier=(49.4274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420436,0.0693541,-5.31309e-05,1.27644e-08,1.87889e-12,45397.5,26.8857], Tmin=(100,'K'), Tmax=(991.762,'K')), NASAPolynomial(coeffs=[17.2933,0.0202381,-7.48502e-06,1.33305e-09,-9.25067e-14,41119.4,-59.0716], Tmin=(991.762,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + other(ketene_1C-C_1C-H) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJC=C=O) + radical(AllylJ2_triplet) + radical(C=COJ)"""),
)

species(
    label = '[CH]C(=C=O)C=[C]C[O](25103)',
    structure = SMILES('[CH]C(=C=O)C=[C]C[O]'),
    E0 = (645.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,350,440,435,1725,3010,987.5,1337.5,450,1655,733.597,733.635,733.655,733.722,733.843,733.911],'cm^-1')),
        HinderedRotor(inertia=(2.42857,'amu*angstrom^2'), symmetry=1, barrier=(55.8375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4283,'amu*angstrom^2'), symmetry=1, barrier=(55.8313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.42821,'amu*angstrom^2'), symmetry=1, barrier=(55.8294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.967786,0.081762,-0.000144716,1.4739e-07,-5.628e-11,77687.1,30.2492], Tmin=(100,'K'), Tmax=(857.712,'K')), NASAPolynomial(coeffs=[0.170225,0.0492398,-2.44586e-05,4.65508e-09,-3.17633e-13,79157,41.7458], Tmin=(857.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + radical(Cds_S) + radical(CCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH][C]=C=O)C[O](25104)',
    structure = SMILES('[CH]=C([CH][C]=C=O)C[O]'),
    E0 = (632.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3025,407.5,1350,352.5,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,2120,512.5,787.5,1685,370,180,180,1446.99],'cm^-1')),
        HinderedRotor(inertia=(0.224352,'amu*angstrom^2'), symmetry=1, barrier=(5.15829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223505,'amu*angstrom^2'), symmetry=1, barrier=(5.13882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.4029,'amu*angstrom^2'), symmetry=1, barrier=(55.2474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594527,0.0878265,-0.000158176,1.51141e-07,-5.46085e-11,76177.5,30.5283], Tmin=(100,'K'), Tmax=(859.562,'K')), NASAPolynomial(coeffs=[5.29204,0.0369041,-1.85964e-05,3.54975e-09,-2.41944e-13,76443.6,14.8216], Tmin=(859.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCCJ=C=O) + radical(CCOJ) + radical(CCJC(C)=C=O)"""),
)

species(
    label = '[CH]C([CH][C]=C=O)=CO(25105)',
    structure = SMILES('[CH]C([CH][C]=C=O)=CO'),
    E0 = (474.686,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,350,440,435,1725,2120,512.5,787.5,1685,370,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.28178,0.0893327,-0.000106395,6.31041e-08,-1.44885e-11,57249.9,29.6235], Tmin=(100,'K'), Tmax=(1076.09,'K')), NASAPolynomial(coeffs=[19.0571,0.0174474,-6.19163e-06,1.02588e-09,-6.64262e-14,53087.8,-65.0935], Tmin=(1076.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + radical(AllylJ2_triplet) + radical(CCJC(C)=C=O) + radical(CCCJ=C=O)"""),
)

species(
    label = '[C-]#[O+](374)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (299.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([180],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33667,0.00896487,-2.66756e-05,3.61071e-08,-1.57199e-11,36069.2,-1.20266], Tmin=(100,'K'), Tmax=(865.594,'K')), NASAPolynomial(coeffs=[-0.394107,0.0117562,-6.47408e-06,1.26375e-09,-8.67562e-14,37256.3,19.3844], Tmin=(865.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.89,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(CsJ2_singlet-CsH)"""),
)

species(
    label = '[CH]=CC=C=C[O](22472)',
    structure = SMILES('[CH]=CC=C=C[O]'),
    E0 = (414.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180],'cm^-1')),
        HinderedRotor(inertia=(1.39997,'amu*angstrom^2'), symmetry=1, barrier=(32.1881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3805.5,'J/mol'), sigma=(6.06656,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=594.41 K, Pc=38.67 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36887,0.0445333,-1.01997e-05,-3.31665e-08,2.01212e-11,49931.8,19.9718], Tmin=(100,'K'), Tmax=(920.489,'K')), NASAPolynomial(coeffs=[19.0433,0.00260849,1.28086e-06,-3.15581e-10,1.86821e-14,45200.3,-71.8594], Tmin=(920.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=COC=C=C[O](22649)',
    structure = SMILES('[CH]=C=COC=C=C[O]'),
    E0 = (338.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,563.333,586.667,610,1970,2140,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.5042,'amu*angstrom^2'), symmetry=1, barrier=(34.5844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52359,'amu*angstrom^2'), symmetry=1, barrier=(35.0304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.00261,0.0837508,-9.25215e-05,4.69213e-08,-8.65362e-12,40868.9,32.6211], Tmin=(100,'K'), Tmax=(1585.25,'K')), NASAPolynomial(coeffs=[24.0008,0.00172769,3.00543e-06,-7.85906e-10,5.8134e-14,35320.5,-92.0229], Tmin=(1585.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C1[CH]OC(C=1)=C[O](25106)',
    structure = SMILES('[CH]C1[CH]OC(C=1)=C[O]'),
    E0 = (331.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,3010,987.5,1337.5,450,1655,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20364,0.0332486,7.59406e-05,-1.39191e-07,6.05897e-11,39941.9,23.6341], Tmin=(100,'K'), Tmax=(924.355,'K')), NASAPolynomial(coeffs=[24.2418,0.00631795,1.56587e-06,-3.90427e-10,1.7823e-14,32574.2,-102.514], Tmin=(924.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(C=CCJ(O)C)"""),
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
    label = '[CH]C(C=O)=CC#C(22634)',
    structure = SMILES('[CH]C(C=O)=CC#C'),
    E0 = (474.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,350,440,435,1725,2175,525,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.0549,'amu*angstrom^2'), symmetry=1, barrier=(47.2462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05334,'amu*angstrom^2'), symmetry=1, barrier=(47.2103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.0545,'amu*angstrom^2'), symmetry=1, barrier=(47.237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (92.0954,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52185,0.0555926,-4.20631e-05,1.56625e-08,-2.41268e-12,57164.2,21.3186], Tmin=(100,'K'), Tmax=(1465.47,'K')), NASAPolynomial(coeffs=[11.8146,0.0274985,-1.33072e-05,2.58108e-09,-1.81074e-13,54147.4,-32.272], Tmin=(1465.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsCtH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=C(C=O)C=C=C[O](25107)',
    structure = SMILES('[C]=C(C=O)C=C=C[O]'),
    E0 = (611.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,540,610,2055,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24244,'amu*angstrom^2'), symmetry=1, barrier=(28.5661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25638,'amu*angstrom^2'), symmetry=1, barrier=(28.8866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0809671,0.078964,-9.09282e-05,4.81121e-08,-9.62937e-12,73680.8,26.3926], Tmin=(100,'K'), Tmax=(1243.18,'K')), NASAPolynomial(coeffs=[22.8535,0.00517025,-1.88898e-06,3.63562e-10,-2.71696e-14,67978.5,-89.2451], Tmin=(1243.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(611.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C1C=C(C=O)C1[O](25108)',
    structure = SMILES('[CH]=C1C=C(C=O)C1[O]'),
    E0 = (388.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924193,0.0547406,-2.54396e-05,-1.16491e-08,9.0735e-12,46852.8,24.7769], Tmin=(100,'K'), Tmax=(1037.76,'K')), NASAPolynomial(coeffs=[18.2795,0.0148297,-6.75473e-06,1.40344e-09,-1.06874e-13,41797.6,-66.5964], Tmin=(1037.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cd-CdCs(CO)) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(C#CC=O)C=O(25109)',
    structure = SMILES('[CH]=C(C#CC=O)C=O'),
    E0 = (273.889,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2695,2870,700,800,1380,1410,450,500,1750,1800,900,1100,350,440,435,1725,2100,2250,500,550],'cm^-1')),
        HinderedRotor(inertia=(1.23805,'amu*angstrom^2'), symmetry=1, barrier=(28.4652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23513,'amu*angstrom^2'), symmetry=1, barrier=(28.3981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24505,'amu*angstrom^2'), symmetry=1, barrier=(28.6262,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (107.087,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13349,0.0628941,-6.22038e-05,2.95365e-08,-5.58207e-12,33044.4,23.6488], Tmin=(100,'K'), Tmax=(1263.5,'K')), NASAPolynomial(coeffs=[14.967,0.0190983,-1.02087e-05,2.10123e-09,-1.5347e-13,29548.8,-46.3252], Tmin=(1263.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtCs) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C=CC=O)C=O(25110)',
    structure = SMILES('[CH]C(=C=CC=O)C=O'),
    E0 = (260.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19669,0.0701325,-7.04506e-05,3.97338e-08,-9.88351e-12,31484.1,23.7005], Tmin=(100,'K'), Tmax=(923.537,'K')), NASAPolynomial(coeffs=[8.15256,0.0400051,-2.15174e-05,4.41045e-09,-3.21419e-13,30199.3,-9.30413], Tmin=(923.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C=O)=CC=C=O(25111)',
    structure = SMILES('[CH]C(C=O)=CC=C=O'),
    E0 = (232.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00836,0.0726713,-9.85105e-05,8.75038e-08,-3.31217e-11,28036.4,26.5863], Tmin=(100,'K'), Tmax=(741.001,'K')), NASAPolynomial(coeffs=[5.18936,0.0423911,-2.16058e-05,4.271e-09,-3.02682e-13,27628.5,9.09709], Tmin=(741.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C=O)C=CC=O(25112)',
    structure = SMILES('[CH]C(=C=O)C=CC=O'),
    E0 = (247.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.933076,0.0757325,-0.000110921,1.03261e-07,-3.94812e-11,29841.7,26.5546], Tmin=(100,'K'), Tmax=(777.256,'K')), NASAPolynomial(coeffs=[4.45899,0.0437307,-2.2421e-05,4.41671e-09,-3.11269e-13,29712.1,13.1251], Tmin=(777.256,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsCs) + group(Cd-Cd(CCO)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1=COC(C=O)=C1(25113)',
    structure = SMILES('[CH]C1=COC(C=O)=C1'),
    E0 = (164.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26009,0.0502016,-1.74187e-05,-8.73782e-09,5.11046e-12,19908.9,23.4377], Tmin=(100,'K'), Tmax=(1183.62,'K')), NASAPolynomial(coeffs=[13.678,0.027506,-1.3078e-05,2.57253e-09,-1.83818e-13,15619.4,-44.2673], Tmin=(1183.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(164.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Furan) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]C=C1[CH]OC=C=C1(25114)',
    structure = SMILES('[O]C=C1[CH]OC=C=C1'),
    E0 = (117.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22044,0.0273263,9.77896e-05,-1.64995e-07,6.95734e-11,14216.5,17.0539], Tmin=(100,'K'), Tmax=(940.714,'K')), NASAPolynomial(coeffs=[27.6409,0.000980083,2.6763e-06,-4.13076e-10,9.6179e-15,5440.61,-129.019], Tmin=(940.714,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + ring(Cyclohexane) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'O=CC1=CC(C=O)=C1(25115)',
    structure = SMILES('O=CC1=CC(C=O)=C1'),
    E0 = (195.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (108.095,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.559466,0.0640819,-5.04012e-05,1.32407e-08,3.19995e-13,23661,20.2127], Tmin=(100,'K'), Tmax=(1134.84,'K')), NASAPolynomial(coeffs=[19.873,0.0136584,-7.08416e-06,1.50005e-09,-1.1301e-13,18140.8,-80.4148], Tmin=(1134.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + ring(cyclobutadiene_13)"""),
)

species(
    label = '[C]=CC(=[CH])C=O(23374)',
    structure = SMILES('[C]=CC(=[CH])C=O'),
    E0 = (785.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,350,440,435,1725,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(1.03003,'amu*angstrom^2'), symmetry=1, barrier=(23.6825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02767,'amu*angstrom^2'), symmetry=1, barrier=(23.6282,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48528,0.0537153,-6.04066e-05,3.24696e-08,-6.79757e-12,94525.6,19.0431], Tmin=(100,'K'), Tmax=(1166.73,'K')), NASAPolynomial(coeffs=[13.9271,0.0110593,-5.56555e-06,1.13322e-09,-8.29042e-14,91622.4,-42.9], Tmin=(1166.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(Cds_P)"""),
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
    E0 = (300.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (390.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (450.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (455.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (559.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (383.246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (404.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (510.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (522.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (492.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (610.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (333.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (372.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (539.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (554.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (615.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (655.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (681.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (487.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (438.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (438.368,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (425.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (364.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (624.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (325.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (452.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (424.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (713.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (591.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (554.965,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (401.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (670.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (657.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (499.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1019.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (651.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (481.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (881.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (823.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (388.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (499.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (487.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (463.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (344.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (356.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (317.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (308.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (852.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['C#CC=O(21959)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[O]C=[C]C1C=C1C=O(25077)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(90.2505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 89.3 to 90.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[O]C=C=CC1=CC1[O](25078)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(9.36651e+10,'s^-1'), n=0.5685, Ea=(150.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonylbond_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic
Ea raised from 149.0 to 150.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]=C1C=C=COC1[O](25079)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.04811e+11,'s^-1'), n=0.222, Ea=(155.084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;multiplebond_intra;radadd_intra_O] + [R7_SMMS;multiplebond_intra;radadd_intra] for rate rule [R7_SMMS;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', '[CH]=C(C=O)C=C=C=O(25080)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C=C[O](8556)', 'C#CC=O(21959)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0669079,'m^3/(mol*s)'), n=2.39465, Ea=(29.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-De_Ct-H;CJ] for rate rule [Ct-CO_Ct-H;CdsJ=Cdd]
Euclidian distance = 1.41421356237
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=O(373)', 'C#CC=C=C[O](23376)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0345031,'m^3/(mol*s)'), n=2.33733, Ea=(26.7584,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cd_Ct-H;CJ] for rate rule [Ct-Cd_Ct-H;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C(C=O)=CC#CO(25081)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH2]C(C#CC=O)=C[O](25082)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['C=C([C]=O)C=C=C[O](25083)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]C(=C=C=CO)C=O(25084)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.63625e+10,'s^-1'), n=1.0925, Ea=(294.328,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R4H_MMS;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH2]C(C=C=C=O)=C[O](25085)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]C(=C=O)C=C=CO(25086)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.33753e+06,'s^-1'), n=1.02312, Ea=(72.6006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H;O_rad_out;XH_out] for rate rule [R6H;O_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C=C[O](8556)', '[CH]=C=C[O](8556)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(81155.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Cd_allenic]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=O(373)', '[CH]=C=C[C]=C[O](23385)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.7839e+06,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;CO_pri_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -3.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH]C(C#CC=O)=C[O](25087)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]C(=C=O)C=C=C[O](25088)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.53107e+07,'m^3/(mol*s)'), n=-0.133333, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_sec_rad;Y_rad] for rate rule [CO_rad/OneDe;H_rad]
Euclidian distance = 1.41421356237
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[CH]C(C=O)=C[C]=C=O(25089)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]=C(C=O)C=C1[CH]O1(25090)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['O=C[C]1[CH]C(C=O)=C1(25091)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[O]C=C=CC1[CH]OC=1(25092)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]=C(C=O)C1[C]=CO1(25093)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[O]C1[C]=CC(C=O)=C1(25094)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.58e+12,'s^-1'), n=-0.292, Ea=(64.2244,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 838 used for R5_MS;doublebond_intra_CdCdd;radadd_intra_cdsingleH
Exact match found for rate rule [R5_MS;doublebond_intra_CdCdd;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]=C1[CH]OOC=C=C1(25095)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(324.35,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 320.3 to 324.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['C=C(C=O)C=C=C=O(25096)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C(=[C]C=C[O])C=O(25097)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C(C=C[C]=O)=C[O](25098)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C(=C=C=C[O])C[O](25099)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C(=C=O)C[C]=C[O](25100)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C(=C=C=C[O])[CH]O(25101)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C(=C=O)[CH]C=C[O](25102)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]C(=C=O)C=[C]C[O](25103)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([CH][C]=C=O)C[O](25104)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C([CH][C]=C=O)=CO(25105)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[C-]#[O+](374)', '[CH]=CC=C=C[O](22472)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C=COC=C=C[O](22649)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C1[CH]OC(C=1)=C[O](25106)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7.69248e+14,'s^-1'), n=-0.917475, Ea=(150.312,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O(T)(63)', '[CH]C(C=O)=CC#C(22634)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(8)', '[C]=C(C=O)C=C=C[O](25107)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]=C1C=C(C=O)C1[O](25108)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(88.2296,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingle] for rate rule [R5_DS_CO;carbonylbond_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.44948974278
family: Intra_R_Add_Exocyclic
Ea raised from 85.4 to 88.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['H(8)', '[CH]=C(C#CC=O)C=O(25109)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2371.94,'m^3/(mol*s)'), n=1.49517, Ea=(13.5032,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-De_Ct-De;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]C(=C=CC=O)C=O(25110)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.7652e+07,'s^-1'), n=1.65613, Ea=(187.464,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_D;Cd_rad_out_single;Cd_H_out_singleDe] + [R2H_D;Cd_rad_out_singleDe;Cd_H_out_single] for rate rule [R2H_D;Cd_rad_out_singleDe;Cd_H_out_singleDe]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]C(C=O)=CC=C=O(25111)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]C(=C=O)C=CC=O(25112)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleDe;CO_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[CH]C1=COC(C=O)=C1(25113)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(8.72e+09,'s^-1'), n=0.186, Ea=(55.7727,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_DS;multiplebond_intra;radadd_intra_cdsingleDe] for rate rule [R5_DS_CO;carbonyl_intra_H;radadd_intra_cdsingleDe]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['[O]C=C1[CH]OC=C=C1(25114)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.59219e+10,'s^-1'), n=0.253963, Ea=(17.5802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=C(C=O)C=C=C[O](22652)'],
    products = ['O=CC1=CC(C=O)=C1(25115)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=O(373)', '[C]=CC(=[CH])C=O(23374)'],
    products = ['[CH]=C(C=O)C=C=C[O](22652)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4955',
    isomers = [
        '[CH]=C(C=O)C=C=C[O](22652)',
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
    label = '4955',
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

