species(
    label = 'C#CC([O])C([CH2])[O](22358)',
    structure = SMILES('C#CC([O])C([CH2])[O]'),
    E0 = (452.682,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3100,440,815,1455,1000,460.595,460.637,460.648],'cm^-1')),
        HinderedRotor(inertia=(0.123984,'amu*angstrom^2'), symmetry=1, barrier=(18.6684,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.123988,'amu*angstrom^2'), symmetry=1, barrier=(18.6688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638613,'amu*angstrom^2'), symmetry=1, barrier=(96.1474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646979,0.0682674,-7.09259e-05,3.39696e-08,-5.29346e-12,54570.7,27.8485], Tmin=(100,'K'), Tmax=(952.482,'K')), NASAPolynomial(coeffs=[16.4129,0.014078,-4.51673e-06,7.37898e-10,-4.87561e-14,51022.1,-50.3078], Tmin=(952.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=C[O](594)',
    structure = SMILES('C=C[O]'),
    E0 = (-25.1807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34719,0.00128739,5.39982e-05,-7.84138e-08,3.24083e-11,-2992.85,8.97297], Tmin=(100,'K'), Tmax=(914.213,'K')), NASAPolynomial(coeffs=[11.726,-0.0014735,2.90737e-06,-5.96989e-10,3.70275e-14,-5941.49,-38.4465], Tmin=(914.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-25.1807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
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
    label = '[CH]=C1OC1C([CH2])[O](24208)',
    structure = SMILES('[CH]=C1OC1C([CH2])[O]'),
    E0 = (474.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674435,0.0558868,-1.39273e-05,-4.07132e-08,2.48221e-11,57236.6,24.5881], Tmin=(100,'K'), Tmax=(927.562,'K')), NASAPolynomial(coeffs=[23.59,0.00156192,1.96752e-06,-4.20157e-10,2.32167e-14,51071.4,-94.5613], Tmin=(927.562,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CJCO) + radical(CC(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OC([CH2])C1[O](24085)',
    structure = SMILES('[CH]=C1OC([CH2])C1[O]'),
    E0 = (423.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.790102,0.0538157,-6.90402e-06,-5.13193e-08,3.05841e-11,51066.7,21.5699], Tmin=(100,'K'), Tmax=(892.715,'K')), NASAPolynomial(coeffs=[23.1494,0.000383793,4.31688e-06,-1.03226e-09,7.22247e-14,45211.6,-94.1979], Tmin=(892.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CC(C)OJ) + radical(Cds_P) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=C1CC([O])C1[O](24209)',
    structure = SMILES('[CH]=C1CC([O])C1[O]'),
    E0 = (462.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.1271,-0.0180329,0.000124643,-1.19872e-07,2.75373e-11,55295.3,-18.2151], Tmin=(100,'K'), Tmax=(1729.92,'K')), NASAPolynomial(coeffs=[77.7674,0.0153692,-6.68923e-05,1.65925e-08,-1.23805e-12,4184.49,-457.326], Tmin=(1729.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_P) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
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
    label = 'C#CC([O])C(=C)[O](22457)',
    structure = SMILES('C#CC([O])C(=C)[O]'),
    E0 = (230.474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,350,440,435,1725,367.566,367.605,367.749],'cm^-1')),
        HinderedRotor(inertia=(0.234357,'amu*angstrom^2'), symmetry=1, barrier=(22.484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234394,'amu*angstrom^2'), symmetry=1, barrier=(22.4842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4234.77,'J/mol'), sigma=(6.72978,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=661.46 K, Pc=31.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707128,0.0626458,-6.76465e-05,3.59801e-08,-7.26644e-12,27846.5,26.7019], Tmin=(100,'K'), Tmax=(1334.73,'K')), NASAPolynomial(coeffs=[16.68,0.00948755,-1.96138e-06,2.02713e-10,-9.04264e-15,24053.8,-53.2049], Tmin=(1334.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CC(=O)C([CH2])[O](24210)',
    structure = SMILES('C#CC(=O)C([CH2])[O]'),
    E0 = (299.163,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,375,552.5,462.5,1710,331.13,332.443],'cm^-1')),
        HinderedRotor(inertia=(0.253492,'amu*angstrom^2'), symmetry=1, barrier=(19.8266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00151896,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00686432,'amu*angstrom^2'), symmetry=1, barrier=(33.5593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42304,0.0581878,-6.65612e-05,4.08543e-08,-1.01147e-11,36072.5,25.9758], Tmin=(100,'K'), Tmax=(979.121,'K')), NASAPolynomial(coeffs=[10.6093,0.0206585,-9.06584e-06,1.70599e-09,-1.18766e-13,34273.7,-18.1487], Tmin=(979.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO) + radical(C=OCOJ)"""),
)

species(
    label = 'CH2(T)(28)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.58,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CC([O])C=O(22474)',
    structure = SMILES('C#CC([O])C=O'),
    E0 = (139.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,474.687],'cm^-1')),
        HinderedRotor(inertia=(0.059994,'amu*angstrom^2'), symmetry=1, barrier=(9.59305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.659976,'amu*angstrom^2'), symmetry=1, barrier=(105.493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87851,0.0410031,-3.92309e-05,1.95951e-08,-3.83768e-12,16840.6,22.268], Tmin=(100,'K'), Tmax=(1295.51,'K')), NASAPolynomial(coeffs=[11.1826,0.0111839,-3.44062e-06,5.26944e-10,-3.24794e-14,14521.5,-24.6743], Tmin=(1295.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
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
    label = 'C#CC([O])C=C(13772)',
    structure = SMILES('C#CC([O])C=C'),
    E0 = (307.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3010,987.5,1337.5,450,1655,217.186,217.21],'cm^-1')),
        HinderedRotor(inertia=(0.885817,'amu*angstrom^2'), symmetry=1, barrier=(29.6435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.885441,'amu*angstrom^2'), symmetry=1, barrier=(29.6432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67511,0.043956,-2.37695e-05,-3.86807e-09,5.69169e-12,37041.1,21.6699], Tmin=(100,'K'), Tmax=(967.656,'K')), NASAPolynomial(coeffs=[12.4618,0.0157303,-5.38113e-06,9.38759e-10,-6.50844e-14,34187.4,-33.9737], Tmin=(967.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
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
    label = '[CH2][CH][O](719)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (361.021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1878.99],'cm^-1')),
        HinderedRotor(inertia=(0.232981,'amu*angstrom^2'), symmetry=1, barrier=(5.35669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03639,0.0272039,-5.17476e-05,5.40082e-08,-2.05139e-11,43449.8,12.3205], Tmin=(100,'K'), Tmax=(879.689,'K')), NASAPolynomial(coeffs=[2.12305,0.0164211,-7.89343e-06,1.47303e-09,-9.88046e-14,44188.4,19.8945], Tmin=(879.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CJCO) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C([O])C=O(1032)',
    structure = SMILES('[CH2]C([O])C=O'),
    E0 = (82.1174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0159141,'amu*angstrom^2'), symmetry=1, barrier=(18.8176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347156,'amu*angstrom^2'), symmetry=1, barrier=(18.8153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37294,0.0300135,-1.22984e-05,-6.17239e-09,4.64176e-12,9940.08,21.4594], Tmin=(100,'K'), Tmax=(1011.88,'K')), NASAPolynomial(coeffs=[9.84794,0.0127268,-4.85036e-06,8.96671e-10,-6.36509e-14,7799.56,-17.7935], Tmin=(1011.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.1174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CJCO)"""),
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
    label = 'C#CC([O])[C]([CH2])O(23482)',
    structure = SMILES('C#CC([O])[C]([CH2])O'),
    E0 = (398.949,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2175,525,3615,1277.5,1000,313.724,314.251],'cm^-1')),
        HinderedRotor(inertia=(0.00170354,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213568,'amu*angstrom^2'), symmetry=1, barrier=(14.865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.210777,'amu*angstrom^2'), symmetry=1, barrier=(14.8585,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16468,'amu*angstrom^2'), symmetry=1, barrier=(81.1901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462635,0.0793052,-0.00010733,7.16916e-08,-1.79023e-11,48108.6,28.1347], Tmin=(100,'K'), Tmax=(803.274,'K')), NASAPolynomial(coeffs=[14.8164,0.0164524,-6.06472e-06,1.01308e-09,-6.4822e-14,45524.3,-39.7016], Tmin=(803.274,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = 'C#CC([O])[C](C)[O](23354)',
    structure = SMILES('C#CC([O])[C](C)[O]'),
    E0 = (417.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,360,370,350,414.23,414.474,417.02],'cm^-1')),
        HinderedRotor(inertia=(0.703553,'amu*angstrom^2'), symmetry=1, barrier=(85.6835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117196,'amu*angstrom^2'), symmetry=1, barrier=(14.3813,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117934,'amu*angstrom^2'), symmetry=1, barrier=(14.382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72166,0.070884,-8.48644e-05,5.25159e-08,-1.27453e-11,50359.3,26.373], Tmin=(100,'K'), Tmax=(1013.01,'K')), NASAPolynomial(coeffs=[14.06,0.0182159,-6.87709e-06,1.19218e-09,-7.9133e-14,47656.9,-38.1493], Tmin=(1013.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=C=C(O)C([CH2])[O](24211)',
    structure = SMILES('[CH]=C=C(O)C([CH2])[O]'),
    E0 = (333.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.196316,0.0776882,-8.95793e-05,4.85243e-08,-9.75134e-12,40228.5,30.3891], Tmin=(100,'K'), Tmax=(1391.31,'K')), NASAPolynomial(coeffs=[21.4749,0.00462022,8.0124e-07,-3.4331e-10,2.86602e-14,35240,-77.5748], Tmin=(1391.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CC(C)OJ) + radical(C=C=CJ) + radical(CJCO)"""),
)

species(
    label = '[CH]=C=C([O])C([CH2])O(23483)',
    structure = SMILES('[CH]=C=C([O])C([CH2])O'),
    E0 = (240.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3615,1277.5,1000,540,610,2055,350,440,435,1725,3000,3100,440,815,1455,1000,1380,1390,370,380,2900,435,327.026],'cm^-1')),
        HinderedRotor(inertia=(0.171131,'amu*angstrom^2'), symmetry=1, barrier=(12.8237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164446,'amu*angstrom^2'), symmetry=1, barrier=(12.8452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173638,'amu*angstrom^2'), symmetry=1, barrier=(12.8259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.592921,0.0707319,-7.82559e-05,3.79901e-08,-5.17126e-12,29059.7,29.2547], Tmin=(100,'K'), Tmax=(882.212,'K')), NASAPolynomial(coeffs=[17.0827,0.0111199,-2.66476e-06,3.37987e-10,-1.89489e-14,25560.5,-51.5747], Tmin=(882.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(CJCO) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C=C([O])C(C)[O](23355)',
    structure = SMILES('[CH]=C=C([O])C(C)[O]'),
    E0 = (259.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,1380,1390,370,380,2900,435,350,440,435,1725,278.841,280.461],'cm^-1')),
        HinderedRotor(inertia=(0.270579,'amu*angstrom^2'), symmetry=1, barrier=(14.8875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267886,'amu*angstrom^2'), symmetry=1, barrier=(14.8844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.494218,0.0664606,-6.99473e-05,3.65874e-08,-7.29359e-12,31326,28.7806], Tmin=(100,'K'), Tmax=(1347.61,'K')), NASAPolynomial(coeffs=[17.2945,0.0112475,-2.53986e-06,2.96959e-10,-1.50908e-14,27283.4,-55.4821], Tmin=(1347.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(CC(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C#CC(O)[C]([CH2])[O](24212)',
    structure = SMILES('C#CC(O)[C]([CH2])[O]'),
    E0 = (398.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462635,0.0793052,-0.00010733,7.16916e-08,-1.79023e-11,48108.6,28.1347], Tmin=(100,'K'), Tmax=(803.274,'K')), NASAPolynomial(coeffs=[14.8164,0.0164524,-6.06472e-06,1.01308e-09,-6.4822e-14,45524.3,-39.7016], Tmin=(803.274,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CJCO) + radical(C2CsJOH)"""),
)

species(
    label = '[C]#CC(O)C([CH2])[O](24213)',
    structure = SMILES('[C]#CC(O)C([CH2])[O]'),
    E0 = (559.464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,387.991,391.051,392.236],'cm^-1')),
        HinderedRotor(inertia=(0.00110264,'amu*angstrom^2'), symmetry=1, barrier=(0.119628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0897722,'amu*angstrom^2'), symmetry=1, barrier=(9.47303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136132,'amu*angstrom^2'), symmetry=1, barrier=(14.9296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.65307,'amu*angstrom^2'), symmetry=1, barrier=(68.6867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479083,0.0761582,-0.00010076,6.78152e-08,-1.759e-11,67415.9,29.3836], Tmin=(100,'K'), Tmax=(1003.03,'K')), NASAPolynomial(coeffs=[14.9282,0.0152312,-4.7027e-06,6.85458e-10,-3.94132e-14,64683.6,-39.54], Tmin=(1003.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C([CH2])O(23484)',
    structure = SMILES('[C]#CC([O])C([CH2])O'),
    E0 = (559.464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3615,1277.5,1000,2175,525,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,387.991,391.051,392.236],'cm^-1')),
        HinderedRotor(inertia=(0.00110264,'amu*angstrom^2'), symmetry=1, barrier=(0.119628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0897722,'amu*angstrom^2'), symmetry=1, barrier=(9.47303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136132,'amu*angstrom^2'), symmetry=1, barrier=(14.9296,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.65307,'amu*angstrom^2'), symmetry=1, barrier=(68.6867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479083,0.0761582,-0.00010076,6.78152e-08,-1.759e-11,67415.9,29.3836], Tmin=(100,'K'), Tmax=(1003.03,'K')), NASAPolynomial(coeffs=[14.9282,0.0152312,-4.7027e-06,6.85458e-10,-3.94132e-14,64683.6,-39.54], Tmin=(1003.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO) + radical(CC(C)OJ) + radical(Acetyl)"""),
)

species(
    label = '[C]#CC([O])C(C)[O](23356)',
    structure = SMILES('[C]#CC([O])C(C)[O]'),
    E0 = (578.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,375.289,375.29,375.294,375.307],'cm^-1')),
        HinderedRotor(inertia=(0.00119683,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162057,'amu*angstrom^2'), symmetry=1, barrier=(16.1983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.711973,'amu*angstrom^2'), symmetry=1, barrier=(71.1618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863833,0.0661125,-7.18905e-05,3.91654e-08,-7.79647e-12,69661.3,27.1799], Tmin=(100,'K'), Tmax=(914.833,'K')), NASAPolynomial(coeffs=[13.8735,0.017527,-5.83263e-06,9.41442e-10,-6.02195e-14,66933.8,-36.3237], Tmin=(914.833,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(578.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Acetyl) + radical(CC(C)OJ) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[C]([CH2])[O](24214)',
    structure = SMILES('C#CC([O])[C]([CH2])[O]'),
    E0 = (629.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3000,3100,440,815,1455,1000,180,633.744,647.281],'cm^-1')),
        HinderedRotor(inertia=(0.739219,'amu*angstrom^2'), symmetry=1, barrier=(16.9961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0453344,'amu*angstrom^2'), symmetry=1, barrier=(13.3378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766845,'amu*angstrom^2'), symmetry=1, barrier=(72.112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.543078,0.0763885,-0.000105747,7.32018e-08,-1.95846e-11,75812.6,27.6474], Tmin=(100,'K'), Tmax=(925.818,'K')), NASAPolynomial(coeffs=[14.8862,0.0144148,-5.33125e-06,8.89581e-10,-5.67772e-14,73156.9,-40.4436], Tmin=(925.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C2CsJOH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])[CH][O](1040)',
    structure = SMILES('[CH2]C([O])[CH][O]'),
    E0 = (395.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,180,1186.82,1189.71],'cm^-1')),
        HinderedRotor(inertia=(0.245926,'amu*angstrom^2'), symmetry=1, barrier=(5.65432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.246963,'amu*angstrom^2'), symmetry=1, barrier=(5.67817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77244,0.0562046,-9.65807e-05,8.93143e-08,-3.15207e-11,47596.6,21.2677], Tmin=(100,'K'), Tmax=(864.007,'K')), NASAPolynomial(coeffs=[5.58235,0.0227355,-1.09916e-05,2.06758e-09,-1.39913e-13,47529.1,6.86351], Tmin=(864.007,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CJCO) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH]=C=C([O])C([CH2])[O](23941)',
    structure = SMILES('[CH]=C=C([O])C([CH2])[O]'),
    E0 = (470.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,350,440,435,1725,359.574,359.593],'cm^-1')),
        HinderedRotor(inertia=(0.166275,'amu*angstrom^2'), symmetry=1, barrier=(15.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166184,'amu*angstrom^2'), symmetry=1, barrier=(15.2521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.48384,0.0700429,-8.43892e-05,4.93105e-08,-1.08995e-11,56771.9,29.4474], Tmin=(100,'K'), Tmax=(1205.74,'K')), NASAPolynomial(coeffs=[17.6449,0.0082411,-1.44528e-06,9.96087e-11,-1.37565e-15,52987.6,-55.0871], Tmin=(1205.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO) + radical(CC(C)OJ) + radical(C=C=CJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[C]#CC([O])C([CH2])[O](23942)',
    structure = SMILES('[C]#CC([O])C([CH2])[O]'),
    E0 = (789.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,3000,3100,440,815,1455,1000,337.182,339.218,339.92,854.498],'cm^-1')),
        HinderedRotor(inertia=(0.00146249,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237619,'amu*angstrom^2'), symmetry=1, barrier=(19.1901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.854982,'amu*angstrom^2'), symmetry=1, barrier=(69.9352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.638434,0.0722047,-9.49926e-05,6.29573e-08,-1.60588e-11,95116.6,28.6198], Tmin=(100,'K'), Tmax=(1008.56,'K')), NASAPolynomial(coeffs=[14.8742,0.0134201,-4.10672e-06,5.95637e-10,-3.42384e-14,92363.3,-39.595], Tmin=(1008.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(789.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(Acetyl) + radical(CJCO) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH2]C([O])C1[C]=CO1(24185)',
    structure = SMILES('[CH2]C([O])C1[C]=CO1'),
    E0 = (446.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849172,0.0478264,1.51379e-05,-7.34967e-08,3.69562e-11,53780,24.6043], Tmin=(100,'K'), Tmax=(928.08,'K')), NASAPolynomial(coeffs=[24.7529,7.42334e-05,2.9838e-06,-5.95604e-10,3.27495e-14,46962.7,-101.757], Tmin=(928.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CJCO) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1OC=[C]C1[O](24124)',
    structure = SMILES('[CH2]C1OC=[C]C1[O]'),
    E0 = (340.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21213,0.0433602,1.66509e-05,-7.18956e-08,3.70833e-11,41125.2,20.1627], Tmin=(100,'K'), Tmax=(894.135,'K')), NASAPolynomial(coeffs=[21.4919,0.00177746,3.97118e-06,-9.75203e-10,6.79382e-14,35534.3,-86.3911], Tmin=(894.135,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.959,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(Cds_S) + radical(CC(C)OJ) + radical(CJC(C)OC)"""),
)

species(
    label = '[O]C1[C]=CCC1[O](24215)',
    structure = SMILES('[O]C1[C]=CCC1[O]'),
    E0 = (394.948,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8177,0.0320555,2.93146e-05,-6.41564e-08,2.73309e-11,47594.2,22.8267], Tmin=(100,'K'), Tmax=(972.283,'K')), NASAPolynomial(coeffs=[15.3546,0.0154848,-5.47439e-06,1.08007e-09,-8.37745e-14,43112.7,-51.6096], Tmin=(972.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.948,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(cyclopentene-vinyl)"""),
)

species(
    label = 'C#CC([O])C(=C)O(23477)',
    structure = SMILES('C#CC([O])C(=C)O'),
    E0 = (92.6688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,350,440,435,1725,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,2175,525,3615,1277.5,1000,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0867,'amu*angstrom^2'), symmetry=1, barrier=(24.9853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0867,'amu*angstrom^2'), symmetry=1, barrier=(24.9854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08744,'amu*angstrom^2'), symmetry=1, barrier=(25.0023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0465299,0.0710705,-7.51927e-05,3.77851e-08,-7.04546e-12,11306.6,27.9137], Tmin=(100,'K'), Tmax=(1520.65,'K')), NASAPolynomial(coeffs=[20.2264,0.00623644,1.11937e-07,-2.05469e-10,1.8465e-14,6471.39,-74.0139], Tmin=(1520.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(92.6688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])C(C)=O(22377)',
    structure = SMILES('C#CC([O])C(C)=O'),
    E0 = (84.4749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,375,552.5,462.5,1710,180,3474.69],'cm^-1')),
        HinderedRotor(inertia=(0.722061,'amu*angstrom^2'), symmetry=1, barrier=(16.6016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0233359,'amu*angstrom^2'), symmetry=1, barrier=(16.5816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14494,'amu*angstrom^2'), symmetry=1, barrier=(82.6416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31117,0.0552658,-5.22137e-05,2.64285e-08,-5.36736e-12,10260.1,25.8952], Tmin=(100,'K'), Tmax=(1191.29,'K')), NASAPolynomial(coeffs=[11.9543,0.0195292,-7.2163e-06,1.24709e-09,-8.28764e-14,7724.25,-27.3148], Tmin=(1191.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.4749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = 'C#CC(O)C(=C)[O](24216)',
    structure = SMILES('C#CC(O)C(=C)[O]'),
    E0 = (0.112853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.585225,0.0661767,-7.20202e-05,3.91439e-08,-8.12145e-12,144.167,27.3301], Tmin=(100,'K'), Tmax=(1285.23,'K')), NASAPolynomial(coeffs=[16.6857,0.0113856,-2.60899e-06,3.04916e-10,-1.52534e-14,-3607.71,-52.8812], Tmin=(1285.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(0.112853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C#CC(=O)C([CH2])O(23478)',
    structure = SMILES('C#CC(=O)C([CH2])O'),
    E0 = (55.4297,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,375,552.5,462.5,1710,2175,525,3615,1277.5,1000,1016.21],'cm^-1')),
        HinderedRotor(inertia=(0.222883,'amu*angstrom^2'), symmetry=1, barrier=(5.12453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.021959,'amu*angstrom^2'), symmetry=1, barrier=(16.0918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.699882,'amu*angstrom^2'), symmetry=1, barrier=(16.0917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.40733,'amu*angstrom^2'), symmetry=1, barrier=(32.3573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15113,0.0615331,-6.44501e-05,3.51398e-08,-7.65546e-12,6770.24,26.1672], Tmin=(100,'K'), Tmax=(1111.94,'K')), NASAPolynomial(coeffs=[12.6861,0.0200379,-8.47272e-06,1.57811e-09,-1.0965e-13,4205.03,-30.7062], Tmin=(1111.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.4297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Ct-CtCs) + group(Ct-CtH) + radical(CJCO)"""),
)

species(
    label = 'C#CC(=O)C(C)[O](23346)',
    structure = SMILES('C#CC(=O)C(C)[O]'),
    E0 = (87.5738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,375,552.5,462.5,1710,283.093,2023.48],'cm^-1')),
        HinderedRotor(inertia=(0.00210447,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.293096,'amu*angstrom^2'), symmetry=1, barrier=(16.6632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.431743,'amu*angstrom^2'), symmetry=1, barrier=(24.5561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56284,0.0532122,-4.78054e-05,2.32074e-08,-4.65143e-12,10620.7,24.8351], Tmin=(100,'K'), Tmax=(1179.25,'K')), NASAPolynomial(coeffs=[10.3046,0.0235596,-1.00867e-05,1.88347e-09,-1.30673e-13,8559.03,-18.7801], Tmin=(1179.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.5738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=OCOJ)"""),
)

species(
    label = 'C#CC([O])C[CH][O](22354)',
    structure = SMILES('C#CC([O])C[CH][O]'),
    E0 = (429.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,3025,407.5,1350,352.5,379.171,379.438,379.734,3275.28],'cm^-1')),
        HinderedRotor(inertia=(0.120234,'amu*angstrom^2'), symmetry=1, barrier=(12.2469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.529125,'amu*angstrom^2'), symmetry=1, barrier=(53.9743,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.530011,'amu*angstrom^2'), symmetry=1, barrier=(53.977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849108,0.0738176,-0.000106349,8.49713e-08,-2.70361e-11,51751.4,27.4389], Tmin=(100,'K'), Tmax=(845.308,'K')), NASAPolynomial(coeffs=[9.70457,0.0261451,-1.15187e-05,2.10869e-09,-1.41931e-13,50460.3,-12.5763], Tmin=(845.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCOJ) + radical(CCsJOH) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])C1CO1(22373)',
    structure = SMILES('C#CC([O])C1CO1'),
    E0 = (200.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.040572,0.0648235,-6.44752e-05,3.27255e-08,-6.07731e-12,24251.7,26.6251], Tmin=(100,'K'), Tmax=(1627.25,'K')), NASAPolynomial(coeffs=[14.3063,0.0121509,3.07095e-07,-4.63836e-10,4.3153e-14,21939.8,-41.9827], Tmin=(1627.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(200.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Ethylene_oxide) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC1OOC1[CH2](22371)',
    structure = SMILES('C#CC1OOC1[CH2]'),
    E0 = (400.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32479,0.0439869,1.24798e-05,-6.74491e-08,3.67268e-11,48337.7,22.6824], Tmin=(100,'K'), Tmax=(860.323,'K')), NASAPolynomial(coeffs=[19.3342,0.00451654,4.12382e-06,-1.17195e-09,8.95931e-14,43600.9,-71.0132], Tmin=(860.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(12dioxetane) + radical(CJCOOH)"""),
)

species(
    label = 'C#CC1OCC1[O](22362)',
    structure = SMILES('C#CC1OCC1[O]'),
    E0 = (194.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66213,0.0373255,2.39879e-05,-7.3441e-08,3.73793e-11,23522.4,20.4886], Tmin=(100,'K'), Tmax=(860.444,'K')), NASAPolynomial(coeffs=[16.7628,0.00844606,2.30037e-06,-8.27296e-10,6.59236e-14,19394.2,-58.982], Tmin=(860.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + ring(Oxetane) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=[C]C([O])C=C(13781)',
    structure = SMILES('[CH]=[C]C([O])C=C'),
    E0 = (626.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,1380,1390,370,380,2900,435,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,373.097,374.84],'cm^-1')),
        HinderedRotor(inertia=(0.143307,'amu*angstrom^2'), symmetry=1, barrier=(14.1956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144824,'amu*angstrom^2'), symmetry=1, barrier=(14.2107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73288,0.0475138,-4.00564e-05,1.72941e-08,-3.02981e-12,75399.7,24.7241], Tmin=(100,'K'), Tmax=(1346.51,'K')), NASAPolynomial(coeffs=[11.3355,0.0189878,-8.27847e-06,1.56051e-09,-1.08621e-13,72813.8,-24.4597], Tmin=(1346.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C#CC([O])[CH][O](23227)',
    structure = SMILES('C#CC([O])[CH][O]'),
    E0 = (453.158,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,750,770,3400,2100,1380,1390,370,380,2900,435,2175,525,270.432,271.639,2343.48],'cm^-1')),
        HinderedRotor(inertia=(0.6648,'amu*angstrom^2'), symmetry=1, barrier=(35.16,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41037,'amu*angstrom^2'), symmetry=1, barrier=(72.9656,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52187,0.0587168,-9.13094e-05,7.51377e-08,-2.39962e-11,54587.5,22.7788], Tmin=(100,'K'), Tmax=(882.675,'K')), NASAPolynomial(coeffs=[8.70407,0.017867,-7.78132e-06,1.39446e-09,-9.1791e-14,53643,-9.14279], Tmin=(882.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = '[CH]=C=CC([CH2])[O](19494)',
    structure = SMILES('[CH]=C=CC([CH2])[O]'),
    E0 = (547.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,540,610,2055,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,188.029],'cm^-1')),
        HinderedRotor(inertia=(0.825126,'amu*angstrom^2'), symmetry=1, barrier=(20.4654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.823216,'amu*angstrom^2'), symmetry=1, barrier=(20.4686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08061,0.0556283,-5.49761e-05,2.75102e-08,-5.32927e-12,65982.7,25.7538], Tmin=(100,'K'), Tmax=(1334.79,'K')), NASAPolynomial(coeffs=[14.7888,0.0122231,-3.58522e-06,5.37588e-10,-3.29651e-14,62530.3,-43.5629], Tmin=(1334.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CJCO) + radical(CC(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]C([O])C([O])C#C(24217)',
    structure = SMILES('[CH]C([O])C([O])C#C'),
    E0 = (689.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2175,525,180,180,623.202,623.227,623.232,623.245],'cm^-1')),
        HinderedRotor(inertia=(0.0736938,'amu*angstrom^2'), symmetry=1, barrier=(20.3127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271579,'amu*angstrom^2'), symmetry=1, barrier=(74.8704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00494632,'amu*angstrom^2'), symmetry=1, barrier=(1.3633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7553,0.0660305,-7.03141e-05,3.39723e-08,-5.28255e-12,83026.2,27.8303], Tmin=(100,'K'), Tmax=(947.095,'K')), NASAPolynomial(coeffs=[16.4527,0.0118534,-3.70457e-06,5.97426e-10,-3.94434e-14,79509.3,-49.9173], Tmin=(947.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(CCJ2_triplet)"""),
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
    E0 = (452.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (537.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (574.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (524.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (463.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (543.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (576.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (550.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (452.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (465.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (692.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (566.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (594.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (606.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (527.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (570.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (527.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (643.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (603.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (618.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (630.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (841.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (947.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (682.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1001.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (694.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (511.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (511.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (475.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (475.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (516.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (516.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (516.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (610,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (458.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (460.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (460.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1033.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (868.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (954.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (901.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C=C[O](594)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH]=C1OC1C([CH2])[O](24208)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_T;triplebond_intra_H;radadd_intra] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH]=C1OC([CH2])C1[O](24085)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_O] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH]=C1CC([O])C1[O](24209)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(177207,'s^-1'), n=1.88643, Ea=(71.6589,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(8)', 'C#CC([O])C(=C)[O](22457)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(8)', 'C#CC(=O)C([CH2])[O](24210)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.59e+07,'cm^3/(mol*s)'), n=1.84, Ea=(32.6352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2830 used for CO-CtCs_O;HJ
Exact match found for rate rule [CO-CtCs_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(28)', 'C#CC([O])C=O(22474)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['O(T)(63)', 'C#CC([O])C=C(13772)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2803 used for Cds-CsH_Cds-HH;O_atom_triplet
Exact match found for rate rule [Cds-CsH_Cds-HH;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[O](594)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(207.988,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 205.0 to 208.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH][O](719)', 'C#CC=O(21959)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.99e+06,'cm^3/(mol*s)'), n=2.12, Ea=(20.0832,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [CO-CtH_O;YJ] for rate rule [CO-CtH_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([O])C=O(1032)', '[C]#C(5143)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CJ] for rate rule [CO-CsH_O;CtJ_Ct]
Euclidian distance = 2.82842712475
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC([O])[C]([CH2])O(23482)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC([O])[C](C)[O](23354)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH]=C=C(O)C([CH2])[O](24211)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.70223e+09,'s^-1'), n=1.15155, Ea=(153.908,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_OneDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_Ct]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH]=C=C([O])C([CH2])O(23483)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH]=C=C([O])C(C)[O](23355)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC(O)[C]([CH2])[O](24212)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[C]#CC(O)C([CH2])[O](24213)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;O_H_out] for rate rule [R4H_TSS;Ct_rad_out;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[C]#CC([O])C([CH2])O(23484)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_TSSS;Ct_rad_out;O_H_out]
Euclidian distance = 2.44948974278
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[C]#CC([O])C(C)[O](23356)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(263079,'s^-1'), n=1.73643, Ea=(39.8993,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_TSSS;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 2.2360679775
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH][O](719)', '[CH]=C=C[O](8556)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(8)', 'C#CC([O])[C]([CH2])[O](24214)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])[CH][O](1040)', '[C]#C(5143)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.34536e+08,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_rad/Ct;Y_rad]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(8)', '[CH]=C=C([O])C([CH2])[O](23941)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[C]#CC([O])C([CH2])[O](23942)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.81e+14,'cm^3/(mol*s)','*|/',3), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 61 used for H_rad;Ct_rad/Ct
Exact match found for rate rule [H_rad;Ct_rad/Ct]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH2]C([O])C1[C]=CO1(24185)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[CH2]C1OC=[C]C1[O](24124)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['[O]C1[C]=CCC1[O](24215)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.47e+11,'s^-1'), n=0.15, Ea=(58.576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_T;triplebond_intra_H;radadd_intra] for rate rule [R5_SS_T;triplebond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC([O])C(=C)O(23477)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC([O])C(C)=O(22377)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC(O)C(=C)[O](24216)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC(=O)C([CH2])O(23478)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC(=O)C(C)[O](23346)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC([O])C[CH][O](22354)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC([O])C1CO1(22373)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC1OOC1[CH2](22371)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.41421356237
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C#CC([O])C([CH2])[O](22358)'],
    products = ['C#CC1OCC1[O](22362)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['O(T)(63)', '[CH]=[C]C([O])C=C(13781)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(T)(28)', 'C#CC([O])[CH][O](23227)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['O(T)(63)', '[CH]=C=CC([CH2])[O](19494)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/OneDeC;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['H(8)', '[CH]C([O])C([O])C#C(24217)'],
    products = ['C#CC([O])C([CH2])[O](22358)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '4649',
    isomers = [
        'C#CC([O])C([CH2])[O](22358)',
    ],
    reactants = [
        ('C=C[O](594)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4649',
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

