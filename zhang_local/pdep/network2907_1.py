species(
    label = '[CH2][C]([O])CC([O])O[O](10222)',
    structure = SMILES('[CH2][C]([O])CC([O])O[O]'),
    E0 = (390.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,180,180,1469.82,1470.27],'cm^-1')),
        HinderedRotor(inertia=(0.265041,'amu*angstrom^2'), symmetry=1, barrier=(6.09381,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265898,'amu*angstrom^2'), symmetry=1, barrier=(6.11353,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265653,'amu*angstrom^2'), symmetry=1, barrier=(6.10788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266347,'amu*angstrom^2'), symmetry=1, barrier=(6.12383,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0883529,0.100085,-0.00018252,1.72045e-07,-6.12442e-11,47149,34.1104], Tmin=(100,'K'), Tmax=(862.061,'K')), NASAPolynomial(coeffs=[6.93095,0.037195,-1.89045e-05,3.61094e-09,-2.4573e-13,47126.3,8.82569], Tmin=(862.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJCO) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C(=C)[O](4273)',
    structure = SMILES('[CH2]C(=C)[O]'),
    E0 = (88.2866,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,510.595],'cm^-1')),
        HinderedRotor(inertia=(0.0480287,'amu*angstrom^2'), symmetry=1, barrier=(8.88265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3365.98,'J/mol'), sigma=(5.64088,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=525.76 K, Pc=42.55 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6374,0.0235792,5.32605e-07,-2.30624e-08,1.26355e-11,10673.5,14.3058], Tmin=(100,'K'), Tmax=(894.06,'K')), NASAPolynomial(coeffs=[10.3562,0.00670937,-7.99446e-07,2.86693e-11,-3.46262e-16,8587.33,-26.0166], Tmin=(894.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(88.2866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
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
    label = 'C=C([O])[CH]C([O])O[O](11086)',
    structure = SMILES('C=C([O])[CH]C([O])O[O]'),
    E0 = (112.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,332.54,332.74,332.819,3076.88],'cm^-1')),
        HinderedRotor(inertia=(0.178066,'amu*angstrom^2'), symmetry=1, barrier=(14.0189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529925,'amu*angstrom^2'), symmetry=1, barrier=(41.582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.528323,'amu*angstrom^2'), symmetry=1, barrier=(41.5829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632398,0.0805422,-0.000112791,8.26687e-08,-2.31121e-11,13681.7,29.189], Tmin=(100,'K'), Tmax=(671.565,'K')), NASAPolynomial(coeffs=[10.8836,0.027804,-1.35792e-05,2.6292e-09,-1.83898e-13,12117.2,-17.5825], Tmin=(671.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(112.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=C(C)OJ) + radical(C=CCJCO) + radical(CCOJ)"""),
)

species(
    label = 'C=C([O])C[C]([O])O[O](4925)',
    structure = SMILES('C=C([O])C[C]([O])O[O]'),
    E0 = (201.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,294.098,294.116,294.149,2447.31],'cm^-1')),
        HinderedRotor(inertia=(0.000813145,'amu*angstrom^2'), symmetry=1, barrier=(9.2325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150403,'amu*angstrom^2'), symmetry=1, barrier=(9.23271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582119,'amu*angstrom^2'), symmetry=1, barrier=(35.7342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.682946,0.0801069,-0.000127288,1.10124e-07,-3.74706e-11,24302.6,32.4742], Tmin=(100,'K'), Tmax=(832.743,'K')), NASAPolynomial(coeffs=[8.94746,0.0288006,-1.39606e-05,2.65835e-09,-1.82458e-13,23328.6,-3.46794], Tmin=(832.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cs_P) + radical(CCOJ)"""),
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
    label = '[CH2][C]([CH2])[O](10271)',
    structure = SMILES('[CH2][C]([CH2])[O]'),
    E0 = (537.173,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,278.503],'cm^-1')),
        HinderedRotor(inertia=(0.00215299,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939305,'amu*angstrom^2'), symmetry=1, barrier=(5.13965,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0694,0.0447257,-6.5608e-05,5.12452e-08,-1.57124e-11,64674.3,16.4544], Tmin=(100,'K'), Tmax=(870.707,'K')), NASAPolynomial(coeffs=[8.27065,0.0130302,-5.47972e-06,9.76923e-10,-6.45299e-14,63716,-11.9064], Tmin=(870.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(CJCO) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = 'C=C([O])C[CH][O](1095)',
    structure = SMILES('C=C([O])C[CH][O]'),
    E0 = (152.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3025,407.5,1350,352.5,420.339,420.351,420.354,2703.21],'cm^-1')),
        HinderedRotor(inertia=(0.0679506,'amu*angstrom^2'), symmetry=1, barrier=(8.51826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0679419,'amu*angstrom^2'), symmetry=1, barrier=(8.51882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4072.67,'J/mol'), sigma=(6.65273,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=636.14 K, Pc=31.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47949,0.0596185,-8.40716e-05,6.86859e-08,-2.26135e-11,18476.4,23.798], Tmin=(100,'K'), Tmax=(831.094,'K')), NASAPolynomial(coeffs=[7.56313,0.0247573,-1.10793e-05,2.0546e-09,-1.39651e-13,17657.9,-3.26687], Tmin=(831.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])[CH]C([O])O[O](11087)',
    structure = SMILES('[CH2]C([O])[CH]C([O])O[O]'),
    E0 = (414.747,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,236.078,236.08,236.081,2193.36],'cm^-1')),
        HinderedRotor(inertia=(3.50415e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153882,'amu*angstrom^2'), symmetry=1, barrier=(6.086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.941436,'amu*angstrom^2'), symmetry=1, barrier=(37.2331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153888,'amu*angstrom^2'), symmetry=1, barrier=(6.08601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.276603,0.0932718,-0.000161014,1.48554e-07,-5.27017e-11,50005.6,36.3154], Tmin=(100,'K'), Tmax=(846.4,'K')), NASAPolynomial(coeffs=[7.46285,0.03599,-1.81698e-05,3.49124e-09,-2.3966e-13,49624.5,7.77888], Tmin=(846.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.747,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCJCOOH) + radical(CJCO) + radical(CC(C)OJ) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2][C]([O])C[C](O)O[O](11088)',
    structure = SMILES('[CH2][C]([O])C[C](O)O[O]'),
    E0 = (370.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.22553,0.103417,-0.000180053,1.58201e-07,-5.30273e-11,44703.1,35.0654], Tmin=(100,'K'), Tmax=(869.344,'K')), NASAPolynomial(coeffs=[11.4369,0.0284469,-1.39295e-05,2.6117e-09,-1.75542e-13,43480.6,-14.9349], Tmin=(869.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CJCO) + radical(Cs_P)"""),
)

species(
    label = '[CH2][C]([O])C[C]([O])OO(3427)',
    structure = SMILES('[CH2][C]([O])C[C]([O])OO'),
    E0 = (444.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,363.333,366.667,370,300,400,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139082,0.104428,-0.000188396,1.75326e-07,-6.21805e-11,53561.2,34.8553], Tmin=(100,'K'), Tmax=(848.947,'K')), NASAPolynomial(coeffs=[8.37481,0.0364526,-1.9063e-05,3.69364e-09,-2.53951e-13,53119.6,1.08791], Tmin=(848.947,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(Cs_P) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C([O])C[C]([O])O[O](3425)',
    structure = SMILES('[CH2]C([O])C[C]([O])O[O]'),
    E0 = (419.578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,237.135,237.137,2280.86,2280.87],'cm^-1')),
        HinderedRotor(inertia=(0.179736,'amu*angstrom^2'), symmetry=1, barrier=(7.17243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00299775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179736,'amu*angstrom^2'), symmetry=1, barrier=(7.17246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.879784,'amu*angstrom^2'), symmetry=1, barrier=(35.1074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294818,0.0919032,-0.00015528,1.41266e-07,-4.97647e-11,50587,34.6332], Tmin=(100,'K'), Tmax=(841.019,'K')), NASAPolynomial(coeffs=[7.97497,0.0350425,-1.76013e-05,3.38336e-09,-2.32673e-13,50014.2,3.18573], Tmin=(841.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(Cs_P) + radical(CC(C)OJ) + radical(ROOJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]([O])[CH]C(O)O[O](11089)',
    structure = SMILES('[CH2][C]([O])[CH]C(O)O[O]'),
    E0 = (365.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.241424,0.104758,-0.000185691,1.65374e-07,-5.59225e-11,44121.7,36.7393], Tmin=(100,'K'), Tmax=(872.138,'K')), NASAPolynomial(coeffs=[10.9115,0.0294178,-1.45118e-05,2.7229e-09,-1.82807e-13,43096.2,-10.2673], Tmin=(872.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(CC(C)OJ) + radical(CCJCOOH) + radical(ROOJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2][C](O)[CH]C([O])O[O](11090)',
    structure = SMILES('[CH2][C](O)[CH]C([O])O[O]'),
    E0 = (361.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0557876,0.106176,-0.000204536,1.96403e-07,-7.00589e-11,43549.8,37.1248], Tmin=(100,'K'), Tmax=(874.86,'K')), NASAPolynomial(coeffs=[6.31929,0.0375688,-1.92482e-05,3.65359e-09,-2.46246e-13,43944.4,15.8514], Tmin=(874.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.014,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(CCJCOOH) + radical(CJCO) + radical(CCOJ) + radical(ROOJ)"""),
)

species(
    label = 'C[C]([O])[CH]C([O])O[O](11091)',
    structure = SMILES('C[C]([O])[CH]C([O])O[O]'),
    E0 = (379.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.270278,0.0968785,-0.000178577,1.7203e-07,-6.23548e-11,45797.7,35.1281], Tmin=(100,'K'), Tmax=(860.013,'K')), NASAPolynomial(coeffs=[5.32321,0.0397584,-2.03142e-05,3.89402e-09,-2.65735e-13,46171.8,18.7408], Tmin=(860.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(379.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH) + radical(CCJCOOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C](O)C[C]([O])O[O](11092)',
    structure = SMILES('[CH2][C](O)C[C]([O])O[O]'),
    E0 = (365.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0412091,0.104852,-0.000198961,1.89322e-07,-6.72061e-11,44131.3,35.4555], Tmin=(100,'K'), Tmax=(873.23,'K')), NASAPolynomial(coeffs=[6.84849,0.0365913,-1.86619e-05,3.54144e-09,-2.389e-13,44327.3,11.1628], Tmin=(873.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(C2CsJOH) + radical(CJCO) + radical(Cs_P) + radical(ROOJ)"""),
)

species(
    label = 'C[C]([O])C[C]([O])O[O](11093)',
    structure = SMILES('C[C]([O])C[C]([O])O[O]'),
    E0 = (384.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.286145,0.0955388,-0.000172949,1.64882e-07,-5.94766e-11,46379.1,33.4543], Tmin=(100,'K'), Tmax=(857.109,'K')), NASAPolynomial(coeffs=[5.84548,0.0387929,-1.9735e-05,3.78359e-09,-2.58533e-13,46557.5,14.0909], Tmin=(857.109,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(Cs_P) + radical(C2CsJOH) + radical(CCOJ) + radical(ROOJ) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][C]([O])[CH]C([O])OO(11094)',
    structure = SMILES('[CH2][C]([O])[CH]C([O])OO'),
    E0 = (439.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1310,387.5,850,1000,360,370,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.155687,0.105777,-0.000194055,1.82511e-07,-6.50718e-11,52979.8,36.5318], Tmin=(100,'K'), Tmax=(852.646,'K')), NASAPolynomial(coeffs=[7.85671,0.0374107,-1.96379e-05,3.80304e-09,-2.61066e-13,52732.2,5.71445], Tmin=(852.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CJCO) + radical(C2CsJOH) + radical(CC(C)OJ) + radical(CCJCOOH) + radical(CCOJ)"""),
)

species(
    label = '[CH2][C]([O])C[CH][O](1105)',
    structure = SMILES('[CH2][C]([O])C[CH][O]'),
    E0 = (547.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,360,370,350,224.73,224.731,2467.76,2467.76],'cm^-1')),
        HinderedRotor(inertia=(0.162139,'amu*angstrom^2'), symmetry=1, barrier=(5.81085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.18136,'amu*angstrom^2'), symmetry=1, barrier=(42.3382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162139,'amu*angstrom^2'), symmetry=1, barrier=(5.81082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.912481,0.0804499,-0.000150196,1.43459e-07,-5.10895e-11,66005.8,26.0227], Tmin=(100,'K'), Tmax=(877.591,'K')), NASAPolynomial(coeffs=[5.37351,0.0307986,-1.522e-05,2.85616e-09,-1.91604e-13,66351.9,11.5156], Tmin=(877.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CJCO) + radical(CCsJOH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2][C]([O])[CH]C([O])O[O](11095)',
    structure = SMILES('[CH2][C]([O])[CH]C([O])O[O]'),
    E0 = (591.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,180,180,1560.36,1560.43],'cm^-1')),
        HinderedRotor(inertia=(0.19757,'amu*angstrom^2'), symmetry=1, barrier=(4.54253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197468,'amu*angstrom^2'), symmetry=1, barrier=(4.54017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197562,'amu*angstrom^2'), symmetry=1, barrier=(4.54233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197556,'amu*angstrom^2'), symmetry=1, barrier=(4.54221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.103765,0.102227,-0.00019882,1.91682e-07,-6.86255e-11,71250.4,36.3599], Tmin=(100,'K'), Tmax=(871.804,'K')), NASAPolynomial(coeffs=[6.23956,0.035801,-1.86769e-05,3.56956e-09,-2.41549e-13,71635.1,15.9418], Tmin=(871.804,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(591.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ) + radical(C2CsJOH) + radical(CCJCOOH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][C]([O])C[C]([O])O[O](11096)',
    structure = SMILES('[CH2][C]([O])C[C]([O])O[O]'),
    E0 = (596.206,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,363.333,366.667,370,300,400,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,180,180,1606.86,1607.36],'cm^-1')),
        HinderedRotor(inertia=(0.220394,'amu*angstrom^2'), symmetry=1, barrier=(5.0673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220303,'amu*angstrom^2'), symmetry=1, barrier=(5.06519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220355,'amu*angstrom^2'), symmetry=1, barrier=(5.0664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220346,'amu*angstrom^2'), symmetry=1, barrier=(5.06619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118283,0.100904,-0.000193251,1.8461e-07,-6.57783e-11,71832,34.6909], Tmin=(100,'K'), Tmax=(869.898,'K')), NASAPolynomial(coeffs=[6.7682,0.0348244,-1.80911e-05,3.45754e-09,-2.34214e-13,72018.2,11.2563], Tmin=(869.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.206,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CJCO) + radical(C2CsJOH) + radical(Cs_P) + radical(CC(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'C=C([O])[CH]C(O)O[O](11097)',
    structure = SMILES('C=C([O])[CH]C(O)O[O]'),
    E0 = (-112.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.207119,0.0843003,-0.000105538,6.69052e-08,-1.66484e-11,-13443.7,29.8391], Tmin=(100,'K'), Tmax=(986.106,'K')), NASAPolynomial(coeffs=[15.8382,0.0208956,-9.09167e-06,1.70281e-09,-1.18315e-13,-16526.6,-45.3534], Tmin=(986.106,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-112.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=CCJCO) + radical(ROOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C[C]([O])OO(3412)',
    structure = SMILES('C=C([O])C[C]([O])OO'),
    E0 = (49.1225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,360,370,350,2950,3100,1380,975,1025,1650,350,440,435,1725,2750,2850,1437.5,1250,1305,750,350,371.336,371.379,371.439,1770.56],'cm^-1')),
        HinderedRotor(inertia=(0.00122193,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0897337,'amu*angstrom^2'), symmetry=1, barrier=(8.78108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0896958,'amu*angstrom^2'), symmetry=1, barrier=(8.7812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.35007,'amu*angstrom^2'), symmetry=1, barrier=(34.2594,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4795.65,'J/mol'), sigma=(7.5684,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=749.07 K, Pc=25.1 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541076,0.0821286,-0.000116475,9.20351e-08,-2.96025e-11,6026.98,32.2332], Tmin=(100,'K'), Tmax=(757.599,'K')), NASAPolynomial(coeffs=[10.2149,0.0310407,-1.53005e-05,2.98409e-09,-2.09802e-13,4561.55,-11.7496], Tmin=(757.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.1225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH2]C([O])CC(=O)O[O](11098)',
    structure = SMILES('[CH2]C([O])CC(=O)O[O]'),
    E0 = (-31.0827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.569333,0.0719301,-7.27057e-05,3.75138e-08,-7.68478e-12,-3611.87,31.0481], Tmin=(100,'K'), Tmax=(1183.55,'K')), NASAPolynomial(coeffs=[15.4188,0.0217435,-9.10006e-06,1.68595e-09,-1.16859e-13,-7126.87,-43.0944], Tmin=(1183.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.0827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-O2s(Cds-O2d)) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + radical(CJCO) + radical(C(=O)OOJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C([O])[CH]C([O])OO(4916)',
    structure = SMILES('C=C([O])[CH]C([O])OO'),
    E0 = (-39.2072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.471742,0.0830004,-0.000104759,7.06021e-08,-1.93007e-11,-4593.13,29.0033], Tmin=(100,'K'), Tmax=(886.206,'K')), NASAPolynomial(coeffs=[12.2489,0.0298442,-1.47887e-05,2.92136e-09,-2.08278e-13,-6680.57,-26.392], Tmin=(886.206,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.2072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJCO)"""),
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
    label = '[CH2][C]([O])CC1OO1(11099)',
    structure = SMILES('[CH2][C]([O])CC1OO1'),
    E0 = (335.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.421715,0.0686036,-7.36294e-05,3.92892e-08,-7.99579e-12,40543.6,28.0815], Tmin=(100,'K'), Tmax=(1314.72,'K')), NASAPolynomial(coeffs=[17.4723,0.0115622,-2.6558e-06,3.11446e-10,-1.57174e-14,36506.7,-57.1452], Tmin=(1314.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(dioxirane) + radical(CJCO) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C1([O])CC([O])O1(2139)',
    structure = SMILES('[CH2]C1([O])CC([O])O1'),
    E0 = (123.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831469,0.0588663,-5.48252e-05,2.74032e-08,-5.2972e-12,14924.6,24.9193], Tmin=(100,'K'), Tmax=(1432.09,'K')), NASAPolynomial(coeffs=[13.4078,0.0170272,-3.97177e-06,4.57265e-10,-2.19311e-14,12010.8,-37.8678], Tmin=(1432.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)OC) + radical(CC(C)(O)OJ) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([CH2])([O])C([O])O[O](10223)',
    structure = SMILES('[CH2]C([CH2])([O])C([O])O[O]'),
    E0 = (407.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,180,246.258,1526.63,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.144677,'amu*angstrom^2'), symmetry=1, barrier=(3.3264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144677,'amu*angstrom^2'), symmetry=1, barrier=(3.3264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144677,'amu*angstrom^2'), symmetry=1, barrier=(3.3264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144677,'amu*angstrom^2'), symmetry=1, barrier=(3.3264,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4813.61,'J/mol'), sigma=(7.82887,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=751.87 K, Pc=22.76 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.383932,0.111591,-0.000207901,1.93512e-07,-6.76937e-11,49135.8,32.7532], Tmin=(100,'K'), Tmax=(869.322,'K')), NASAPolynomial(coeffs=[8.89764,0.0352413,-1.81124e-05,3.44945e-09,-2.33329e-13,48793.3,-3.41372], Tmin=(869.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CJC(C)2O) + radical(ROOJ) + radical(CJC(C)2O)"""),
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
    label = '[CH2]C([O])=C[CH][O](4707)',
    structure = SMILES('[CH2]C([O])=C[CH][O]'),
    E0 = (242.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,440,435,1725,548.079,548.792,551.945],'cm^-1')),
        HinderedRotor(inertia=(0.000558036,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000551953,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85419,0.0417783,-2.96654e-05,6.5664e-09,1.15183e-12,29203,24.1344], Tmin=(100,'K'), Tmax=(1002.65,'K')), NASAPolynomial(coeffs=[11.1938,0.0150611,-5.46764e-06,9.6426e-10,-6.61676e-14,26800.2,-23.5914], Tmin=(1002.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=CCJO) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][C]([O])CC1OOO1(11100)',
    structure = SMILES('[CH2][C]([O])CC1OOO1'),
    E0 = (377.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647859,0.0607959,-3.14355e-05,-8.40408e-09,8.49862e-12,45494.8,30.3694], Tmin=(100,'K'), Tmax=(1026.69,'K')), NASAPolynomial(coeffs=[18.648,0.0175163,-7.43019e-06,1.47934e-09,-1.10159e-13,40383.6,-63.8364], Tmin=(1026.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1([O])CC(O[O])O1(10226)',
    structure = SMILES('[CH2]C1([O])CC(O[O])O1'),
    E0 = (120.869,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0128861,0.0755769,-8.10831e-05,4.42697e-08,-9.13154e-12,14691.2,28.511], Tmin=(100,'K'), Tmax=(1360.3,'K')), NASAPolynomial(coeffs=[17.3774,0.0142841,-2.21277e-06,8.67257e-11,4.76147e-15,10913.7,-57.126], Tmin=(1360.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CJC(C)OC) + radical(ROOJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = '[CH2]C1([O])CC([O])OO1(11101)',
    structure = SMILES('[CH2]C1([O])CC([O])OO1'),
    E0 = (98.4162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03788,0.0525348,-8.17534e-06,-3.76504e-08,2.2081e-11,11955.2,27.2115], Tmin=(100,'K'), Tmax=(903.552,'K')), NASAPolynomial(coeffs=[17.2214,0.0151045,-2.83557e-06,3.17623e-10,-1.9457e-14,7634.06,-56.9519], Tmin=(903.552,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(98.4162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(CJCOOH) + radical(CCOJ) + radical(CC(C)(O)OJ)"""),
)

species(
    label = 'C=C([O])CC([O])O[O](4917)',
    structure = SMILES('C=C([O])CC([O])O[O]'),
    E0 = (-4.11921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (117.08,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.679966,0.0789492,-0.000115274,9.57779e-08,-3.21376e-11,-381.544,31.7985], Tmin=(100,'K'), Tmax=(800.282,'K')), NASAPolynomial(coeffs=[9.01088,0.031348,-1.48793e-05,2.83725e-09,-1.96126e-13,-1524.06,-5.34457], Tmin=(800.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-4.11921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH2][C]([O])CC([O])[O](2857)',
    structure = SMILES('[CH2][C]([O])CC([O])[O]'),
    E0 = (393.155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,360,370,350,180,180,1659.15,1659.63,1660.19],'cm^-1')),
        HinderedRotor(inertia=(0.192533,'amu*angstrom^2'), symmetry=1, barrier=(4.42671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192046,'amu*angstrom^2'), symmetry=1, barrier=(4.41551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192844,'amu*angstrom^2'), symmetry=1, barrier=(4.43386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865427,0.0839802,-0.000158955,1.59577e-07,-5.97249e-11,47384,29.9669], Tmin=(100,'K'), Tmax=(856.239,'K')), NASAPolynomial(coeffs=[2.64217,0.040437,-2.09332e-05,4.0421e-09,-2.77259e-13,48371.7,29.2151], Tmin=(856.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(393.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCOJ) + radical(CJCO) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2][C][O](2305)',
    structure = SMILES('[CH2][C][O]'),
    E0 = (641.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1736.9],'cm^-1')),
        HinderedRotor(inertia=(0.485156,'amu*angstrom^2'), symmetry=1, barrier=(11.1547,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.0699,0.0246859,-4.66718e-05,4.62565e-08,-1.71212e-11,77209.5,11.3828], Tmin=(100,'K'), Tmax=(857.616,'K')), NASAPolynomial(coeffs=[3.88751,0.0110738,-5.72565e-06,1.10472e-09,-7.56931e-14,77429.6,9.66472], Tmin=(857.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CH2_triplet) + radical(CCOJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])O[O](607)',
    structure = SMILES('[CH2]C([O])O[O]'),
    E0 = (206.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1967.16],'cm^-1')),
        HinderedRotor(inertia=(0.24983,'amu*angstrom^2'), symmetry=1, barrier=(5.74408,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249705,'amu*angstrom^2'), symmetry=1, barrier=(5.74121,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (75.0434,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29491,0.0444229,-7.76836e-05,7.63138e-08,-2.8547e-11,24875.9,21.4728], Tmin=(100,'K'), Tmax=(845.116,'K')), NASAPolynomial(coeffs=[3.56685,0.022769,-1.15016e-05,2.21713e-09,-1.52683e-13,25219.2,18.8536], Tmin=(845.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCOJ) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][C]([O])C[CH]O[O](11102)',
    structure = SMILES('[CH2][C]([O])C[CH]O[O]'),
    E0 = (554.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,492.5,1135,1000,360,370,350,3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,221.439,221.44,2350.47],'cm^-1')),
        HinderedRotor(inertia=(0.159411,'amu*angstrom^2'), symmetry=1, barrier=(5.54699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65081,'amu*angstrom^2'), symmetry=1, barrier=(57.4428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159413,'amu*angstrom^2'), symmetry=1, barrier=(5.54699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159412,'amu*angstrom^2'), symmetry=1, barrier=(5.54699,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.364669,0.0910888,-0.000161364,1.4639e-07,-5.02233e-11,66759.3,29.2493], Tmin=(100,'K'), Tmax=(878.86,'K')), NASAPolynomial(coeffs=[8.3593,0.0298267,-1.43474e-05,2.66336e-09,-1.77706e-13,66314.8,-2.8223], Tmin=(878.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(CJCO) + radical(CC(C)OJ) + radical(CCsJOOH) + radical(ROOJ)"""),
)

species(
    label = '[CH][C]([O])CC([O])O[O](11103)',
    structure = SMILES('[CH][C]([O])CC([O])O[O]'),
    E0 = (627.586,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (116.072,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.187916,0.097955,-0.000182293,1.72559e-07,-6.1454e-11,75604.9,34.1233], Tmin=(100,'K'), Tmax=(862.956,'K')), NASAPolynomial(coeffs=[7.00703,0.0349078,-1.80558e-05,3.46177e-09,-2.35691e-13,75598.6,9.01274], Tmin=(862.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.586,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCJ2_triplet) + radical(CCOJ) + radical(ROOJ) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
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
    E0 = (390.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (390.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (459.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (390.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (398.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (390.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (551.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (548.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (579.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (535.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (466.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (546.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (535.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (514.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (522.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (523.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (539.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (766.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (803.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (808.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (454.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (454.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (454.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (415.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (390.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (580.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (473.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (564.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (530.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (399.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (399.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (398.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (390.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (641.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (882.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (960.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (839.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2]C(=C)[O](4273)', '[O]OC=O(5472)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', 'C=C([O])[CH]C([O])O[O](11086)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(66.3574,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 62.3 to 66.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=C([O])C[C]([O])O[O](4925)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(=C)[O](4273)', '[O][CH]O[O](8201)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.246938,'m^3/(mol*s)'), n=2.00579, Ea=(73.7851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 70.3 to 73.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]OC=O(5472)', '[CH2][C]([CH2])[O](10271)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(13.7458,'m^3/(mol*s)'), n=1.39198, Ea=(57.1551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;CJ] + [CO-NdH_O;YJ] for rate rule [CO-NdH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O2(2)', 'C=C([O])C[CH][O](1095)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(0.0849,'cm^3/(mol*s)'), n=3.486, Ea=(246.686,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [CO-CsH_O;OJ] for rate rule [CO-CsH_O;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 244.2 to 246.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C([O])[CH]C([O])O[O](11087)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2][C]([O])C[C](O)O[O](11088)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]([O])C[C]([O])OO(3427)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([O])C[C]([O])O[O](3425)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2][C]([O])[CH]C(O)O[O](11089)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2][C](O)[CH]C([O])O[O](11090)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.30814e+09,'s^-1'), n=1.19923, Ea=(155.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['C[C]([O])[CH]C([O])O[O](11091)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2][C](O)C[C]([O])O[O](11092)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(520772,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['C[C]([O])C[C]([O])O[O](11093)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]([O])[CH]C([O])OO(11094)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[CH2][C]([O])C[CH][O](1105)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.18266e+06,'m^3/(mol*s)'), n=0.193158, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;O2_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -25.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][CH]O[O](8201)', '[CH2][C]([CH2])[O](10271)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(8)', '[CH2][C]([O])[CH]C([O])O[O](11095)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(8)', '[CH2][C]([O])C[C]([O])O[O](11096)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['C=C([O])[CH]C(O)O[O](11097)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['C=C([O])C[C]([O])OO(3412)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2]C([O])CC(=O)O[O](11098)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['C=C([O])[CH]C([O])OO(4916)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['C=C([O])C[CH][O](1095)', 'O2(S)(5486)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['O(T)(63)', '[CH2][C]([O])CC1OO1(11099)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(189.53,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['O(T)(63)', '[CH2]C1([O])CC([O])O1(2139)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;Y_rad_intra;OO] for rate rule [R3OO_SS;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([CH2])([O])C([O])O[O](10223)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[O]O(16)', '[CH2]C([O])=C[CH][O](4707)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.91635e+11,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2][C]([O])CC1OOO1(11100)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2]C1([O])CC(O[O])O1(10226)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['[CH2]C1([O])CC([O])OO1(11101)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][C]([O])CC([O])O[O](10222)'],
    products = ['C=C([O])CC([O])O[O](4917)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O(T)(63)', '[CH2][C]([O])CC([O])[O](2857)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][C][O](2305)', '[CH2]C([O])O[O](607)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['O(T)(63)', '[CH2][C]([O])C[CH]O[O](11102)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(8)', '[CH][C]([O])CC([O])O[O](11103)'],
    products = ['[CH2][C]([O])CC([O])O[O](10222)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '2907',
    isomers = [
        '[CH2][C]([O])CC([O])O[O](10222)',
    ],
    reactants = [
        ('[CH2]C(=C)[O](4273)', '[O]OC=O(5472)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '2907',
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

