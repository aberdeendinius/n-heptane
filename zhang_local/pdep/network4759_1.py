species(
    label = 'C=[C]OC=C=C[O](22462)',
    structure = SMILES('C=[C]OC=C=C[O]'),
    E0 = (282.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25544,'amu*angstrom^2'), symmetry=1, barrier=(28.8649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25328,'amu*angstrom^2'), symmetry=1, barrier=(28.8154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.98208,0.0555322,-3.43713e-05,-7.90864e-09,1.06824e-11,34132.4,27.1296], Tmin=(100,'K'), Tmax=(931.667,'K')), NASAPolynomial(coeffs=[18.6352,0.0073376,-1.20845e-06,1.5465e-10,-1.26314e-14,29645.4,-63.2146], Tmin=(931.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO)"""),
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
    label = 'C=C1OC1[C]=C[O](23793)',
    structure = SMILES('C=C1OC1[C]=C[O]'),
    E0 = (252.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1185,0.0414355,2.40097e-05,-8.15926e-08,3.99709e-11,30520,21.1528], Tmin=(100,'K'), Tmax=(923.331,'K')), NASAPolynomial(coeffs=[24.6579,-0.00362618,4.75453e-06,-9.31534e-10,5.58132e-14,23746.9,-103.671], Tmin=(923.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Cds_S) + radical(C=COJ)"""),
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
    label = 'C=[C]OC=C=C=O(23794)',
    structure = SMILES('C=[C]OC=C=C=O'),
    E0 = (326.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,3010,987.5,1337.5,450,1655,2120,512.5,787.5,1685,370,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.248198,'amu*angstrom^2'), symmetry=1, barrier=(5.70657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249856,'amu*angstrom^2'), symmetry=1, barrier=(5.74469,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25512,0.0320431,-1.6297e-05,-5.58136e-09,5.46026e-12,39385.9,9.76221], Tmin=(100,'K'), Tmax=(974.119,'K')), NASAPolynomial(coeffs=[11.4761,0.00872709,-2.7952e-06,5.09322e-10,-3.74582e-14,36899.2,-38.0245], Tmin=(974.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = 'C#COC=C=C[O](23795)',
    structure = SMILES('C#COC=C=C[O]'),
    E0 = (252.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,2995,3025,975,1000,1300,1375,400,500,1630,1680,2175,525,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.53098,'amu*angstrom^2'), symmetry=1, barrier=(35.2003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52577,'amu*angstrom^2'), symmetry=1, barrier=(35.0805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.045816,0.0688729,-7.5888e-05,3.7924e-08,-6.93703e-12,30518.3,25.1592], Tmin=(100,'K'), Tmax=(1567.29,'K')), NASAPolynomial(coeffs=[21.5893,0.000623621,1.90415e-06,-4.71425e-10,3.37549e-14,25337.3,-83.8331], Tmin=(1567.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ)"""),
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
    label = 'C=[C]OC=C=[C]O(23796)',
    structure = SMILES('C=[C]OC=C=[C]O'),
    E0 = (381.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.96132,'amu*angstrom^2'), symmetry=1, barrier=(22.1026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961623,'amu*angstrom^2'), symmetry=1, barrier=(22.1096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96096,'amu*angstrom^2'), symmetry=1, barrier=(22.0944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.431612,0.0684098,-7.93398e-05,4.39539e-08,-9.14313e-12,45971.5,30.8208], Tmin=(100,'K'), Tmax=(1306.49,'K')), NASAPolynomial(coeffs=[18.5685,0.00635644,-6.04321e-07,-4.51975e-11,7.63835e-15,41789.2,-59.3969], Tmin=(1306.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC=C=C[O](23797)',
    structure = SMILES('[CH]=COC=C=C[O]'),
    E0 = (290.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.37214,'amu*angstrom^2'), symmetry=1, barrier=(31.5481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.37227,'amu*angstrom^2'), symmetry=1, barrier=(31.5512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700006,0.0572671,-2.39735e-05,-2.96621e-08,2.1209e-11,35031.1,25.4173], Tmin=(100,'K'), Tmax=(915.612,'K')), NASAPolynomial(coeffs=[23.0289,0.000357221,2.68463e-06,-5.98645e-10,3.81898e-14,29238.8,-89.6402], Tmin=(915.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CO[C]=C=C[O](23798)',
    structure = SMILES('C=CO[C]=C=C[O]'),
    E0 = (282.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.25544,'amu*angstrom^2'), symmetry=1, barrier=(28.8649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25328,'amu*angstrom^2'), symmetry=1, barrier=(28.8154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.98208,0.0555322,-3.43713e-05,-7.90864e-09,1.06824e-11,34132.4,27.1296], Tmin=(100,'K'), Tmax=(931.667,'K')), NASAPolynomial(coeffs=[18.6352,0.0073376,-1.20845e-06,1.5465e-10,-1.26314e-14,29645.4,-63.2146], Tmin=(931.667,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]O[C]=C=CO(23799)',
    structure = SMILES('C=[C]O[C]=C=CO'),
    E0 = (381.088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2950,3100,1380,975,1025,1650,540,610,2055,3615,1277.5,1000,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.96132,'amu*angstrom^2'), symmetry=1, barrier=(22.1026,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.961623,'amu*angstrom^2'), symmetry=1, barrier=(22.1096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96096,'amu*angstrom^2'), symmetry=1, barrier=(22.0944,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.431612,0.0684098,-7.93398e-05,4.39539e-08,-9.14313e-12,45971.5,30.8208], Tmin=(100,'K'), Tmax=(1306.49,'K')), NASAPolynomial(coeffs=[18.5685,0.00635644,-6.04321e-07,-4.51975e-11,7.63835e-15,41789.2,-59.3969], Tmin=(1306.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=CO[CH][C]=C=O(23800)',
    structure = SMILES('C=CO[CH][C]=C=O'),
    E0 = (182.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700545,0.0686723,-8.08194e-05,4.71297e-08,-1.0626e-11,22050.1,24.0647], Tmin=(100,'K'), Tmax=(1095.49,'K')), NASAPolynomial(coeffs=[15.9281,0.0130715,-4.68826e-06,7.99812e-10,-5.31146e-14,18713.7,-50.7883], Tmin=(1095.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C#CO[CH][C]=CO(23801)',
    structure = SMILES('C#CO[CH][C]=CO'),
    E0 = (273.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,2175,525,3615,1277.5,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.33335,'amu*angstrom^2'), symmetry=1, barrier=(30.6562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33341,'amu*angstrom^2'), symmetry=1, barrier=(30.6578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33495,'amu*angstrom^2'), symmetry=1, barrier=(30.693,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.33265,'amu*angstrom^2'), symmetry=1, barrier=(30.6401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.443628,0.0774222,-8.83445e-05,4.56786e-08,-8.63898e-12,33055.9,26.4662], Tmin=(100,'K'), Tmax=(1512.34,'K')), NASAPolynomial(coeffs=[23.4394,6.27192e-06,2.5709e-06,-6.27793e-10,4.55202e-14,27461.4,-93.2481], Tmin=(1512.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
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
    label = 'C=[C]O[C]=C=C[O](23802)',
    structure = SMILES('C=[C]O[C]=C=C[O]'),
    E0 = (522.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03613,'amu*angstrom^2'), symmetry=1, barrier=(23.8228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03814,'amu*angstrom^2'), symmetry=1, barrier=(23.8689,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09568,0.0584211,-6.51691e-05,3.58571e-08,-7.60854e-12,62957.6,29.8503], Tmin=(100,'K'), Tmax=(1164.8,'K')), NASAPolynomial(coeffs=[15.0303,0.0105686,-3.54558e-06,5.87026e-10,-3.85238e-14,59711.4,-39.5018], Tmin=(1164.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=CJO) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]O[CH][C]=C=O(23803)',
    structure = SMILES('C=[C]O[CH][C]=C=O'),
    E0 = (422.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,2120,512.5,787.5,2950,3100,1380,975,1025,1650,349.015,349.028,349.378,349.857],'cm^-1')),
        HinderedRotor(inertia=(0.24624,'amu*angstrom^2'), symmetry=1, barrier=(21.3224,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144108,'amu*angstrom^2'), symmetry=1, barrier=(12.5302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777583,'amu*angstrom^2'), symmetry=1, barrier=(67.5692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00679,0.0692662,-0.000103551,8.05251e-08,-2.45961e-11,50867,26.0966], Tmin=(100,'K'), Tmax=(839.546,'K')), NASAPolynomial(coeffs=[11.2751,0.018078,-8.0468e-06,1.47307e-09,-9.89635e-14,49222.7,-21.1705], Tmin=(839.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C#CO[CH][C]=C[O](23804)',
    structure = SMILES('C#CO[CH][C]=C[O]'),
    E0 = (414.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,750,770,3400,2100,3010,987.5,1337.5,450,1655,2175,525,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.54938,'amu*angstrom^2'), symmetry=1, barrier=(35.6232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.54803,'amu*angstrom^2'), symmetry=1, barrier=(35.5924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5475,'amu*angstrom^2'), symmetry=1, barrier=(35.5802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.223953,0.0673943,-7.40337e-05,3.73788e-08,-7.00752e-12,50041.8,25.4827], Tmin=(100,'K'), Tmax=(1479.82,'K')), NASAPolynomial(coeffs=[20.8922,0.00274445,4.0036e-07,-1.64844e-10,1.26417e-14,44886.4,-79.0806], Tmin=(1479.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]OC=C1[CH]O1(23805)',
    structure = SMILES('C=[C]OC=C1[CH]O1'),
    E0 = (308.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36796,0.0347418,4.15334e-05,-9.8809e-08,4.60505e-11,37279.1,24.3264], Tmin=(100,'K'), Tmax=(919.3,'K')), NASAPolynomial(coeffs=[23.9222,-0.00313107,4.99906e-06,-1.00647e-09,6.15056e-14,30585.8,-96.4369], Tmin=(919.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=CJO) + radical(C=CCJO)"""),
)

species(
    label = 'C=C1O[CH]C1=C[O](23806)',
    structure = SMILES('C=C1O[CH]C1=C[O]'),
    E0 = (86.5843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49771,0.0364367,2.69475e-05,-7.27534e-08,3.35679e-11,10520.9,19.183], Tmin=(100,'K'), Tmax=(939.398,'K')), NASAPolynomial(coeffs=[19.5701,0.00656079,-5.19213e-07,8.63828e-11,-1.40085e-14,5048.2,-77.9323], Tmin=(939.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.5843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=[C]OC1[C]=CO1(23807)',
    structure = SMILES('C=[C]OC1[C]=CO1'),
    E0 = (346.547,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975734,0.0440451,2.11243e-05,-8.08895e-08,4.01722e-11,41809.7,22.7835], Tmin=(100,'K'), Tmax=(923.668,'K')), NASAPolynomial(coeffs=[25.5975,-0.00431342,5.03261e-06,-9.79326e-10,5.87701e-14,34775.7,-107.502], Tmin=(923.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OC=[C]C1[O](23808)',
    structure = SMILES('C=C1OC=[C]C1[O]'),
    E0 = (286.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28148,0.0169326,6.97181e-05,-1.09667e-07,4.49724e-11,34490.4,22.5631], Tmin=(100,'K'), Tmax=(942.372,'K')), NASAPolynomial(coeffs=[16.764,0.00848775,-1.24586e-06,2.46702e-10,-2.72967e-14,29406.2,-58.9401], Tmin=(942.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(Cds_S) + radical(CC(C)OJ)"""),
)

species(
    label = 'C=COC=C=C=O(23809)',
    structure = SMILES('C=COC=C=C=O'),
    E0 = (87.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05071,0.0302223,1.07987e-05,-4.46322e-08,2.17912e-11,10564.6,7.36722], Tmin=(100,'K'), Tmax=(934.884,'K')), NASAPolynomial(coeffs=[15.5494,0.00471344,-1.23821e-08,-2.73621e-11,-2.97426e-15,6631.48,-64.3842], Tmin=(934.884,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=[C]O[C]C=C[O](23810)',
    structure = SMILES('C=[C]O[C]C=C[O]'),
    E0 = (557.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03107,'amu*angstrom^2'), symmetry=1, barrier=(23.7063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03323,'amu*angstrom^2'), symmetry=1, barrier=(23.756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03263,'amu*angstrom^2'), symmetry=1, barrier=(23.7422,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.436024,0.0726473,-7.82889e-05,3.80969e-08,-6.74862e-12,67227.5,30.8615], Tmin=(100,'K'), Tmax=(1640.67,'K')), NASAPolynomial(coeffs=[22.7243,0.000321968,2.3353e-06,-5.55584e-10,3.89092e-14,61762.4,-85.8352], Tmin=(1640.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(CH2_triplet) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]O[CH]C=[C][O](23811)',
    structure = SMILES('C=[C]O[CH]C=[C][O]'),
    E0 = (447.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,221.666,560.445,560.464,560.478,560.573],'cm^-1')),
        HinderedRotor(inertia=(0.0750038,'amu*angstrom^2'), symmetry=1, barrier=(16.7207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749866,'amu*angstrom^2'), symmetry=1, barrier=(16.7208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188404,'amu*angstrom^2'), symmetry=1, barrier=(41.9975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19814,0.05464,-4.63166e-05,1.58932e-08,-8.37517e-13,53886.6,30.3053], Tmin=(100,'K'), Tmax=(1021.91,'K')), NASAPolynomial(coeffs=[14.597,0.0144849,-5.41668e-06,9.8078e-10,-6.86608e-14,50506.4,-37.7668], Tmin=(1021.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(C=CJO) + radical(C=CJO) + radical(C=COJ)"""),
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
    label = '[O]C=C=C[O](22349)',
    structure = SMILES('[O]C=C=C[O]'),
    E0 = (48.0679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10387,0.0254095,2.29542e-05,-6.61705e-08,3.24893e-11,5864.8,15.5654], Tmin=(100,'K'), Tmax=(907.068,'K')), NASAPolynomial(coeffs=[19.4725,-0.00784929,6.29342e-06,-1.25746e-09,8.23938e-14,931.198,-76.3609], Tmin=(907.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.0679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = '[CH]=C=CO[C]=C(19367)',
    structure = SMILES('[CH]=C=CO[C]=C'),
    E0 = (504.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,3010,987.5,1337.5,450,1655,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21977,'amu*angstrom^2'), symmetry=1, barrier=(28.0448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21614,'amu*angstrom^2'), symmetry=1, barrier=(27.9614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34459,0.0510081,-5.26791e-05,2.77943e-08,-5.63687e-12,60792.9,25.1839], Tmin=(100,'K'), Tmax=(1318.63,'K')), NASAPolynomial(coeffs=[13.3316,0.0108536,-2.68749e-06,3.38623e-10,-1.80043e-14,57961.3,-34.7116], Tmin=(1318.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]OC#CC=O(23812)',
    structure = SMILES('C=[C]OC#CC=O'),
    E0 = (288.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2782.5,750,1395,475,1775,1000,2100,2250,500,550,2950,3100,1380,975,1025,1650,322.132,322.169,322.188],'cm^-1')),
        HinderedRotor(inertia=(0.257746,'amu*angstrom^2'), symmetry=1, barrier=(18.9825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.468966,'amu*angstrom^2'), symmetry=1, barrier=(34.5466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257746,'amu*angstrom^2'), symmetry=1, barrier=(18.9806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83893,0.0510805,-5.54803e-05,3.29423e-08,-8.14915e-12,34824.4,24.2214], Tmin=(100,'K'), Tmax=(961.381,'K')), NASAPolynomial(coeffs=[8.6648,0.0226803,-1.11687e-05,2.21462e-09,-1.58661e-13,33511.9,-8.44063], Tmin=(961.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtOs) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]O[C]=CC=O(23813)',
    structure = SMILES('C=[C]O[C]=CC=O'),
    E0 = (325.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,259,259.966,262.465],'cm^-1')),
        HinderedRotor(inertia=(0.00249358,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269305,'amu*angstrom^2'), symmetry=1, barrier=(12.5073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264924,'amu*angstrom^2'), symmetry=1, barrier=(12.5089,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22768,0.0478849,-2.40015e-05,-5.77903e-08,7.41033e-11,39240.9,26.8948], Tmin=(100,'K'), Tmax=(450.809,'K')), NASAPolynomial(coeffs=[5.85399,0.0305626,-1.57876e-05,3.15166e-09,-2.24728e-13,38763,10.615], Tmin=(450.809,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]O[CH]C=C=O(23814)',
    structure = SMILES('C=[C]O[CH]C=C=O'),
    E0 = (184.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11569,0.0637308,-7.52657e-05,4.63108e-08,-1.12963e-11,22260.4,25.2346], Tmin=(100,'K'), Tmax=(1001.81,'K')), NASAPolynomial(coeffs=[12.3943,0.0186979,-7.83828e-06,1.44036e-09,-9.89837e-14,20000.6,-29.1984], Tmin=(1001.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C#CO[CH]C=C[O](23815)',
    structure = SMILES('C#CO[CH]C=C[O]'),
    E0 = (176.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813399,0.0563642,-2.73394e-05,-1.95312e-08,1.54409e-11,21414,22.8851], Tmin=(100,'K'), Tmax=(941.422,'K')), NASAPolynomial(coeffs=[20.8596,0.0053757,-5.66975e-07,8.2277e-11,-1.06878e-14,16124.7,-80.662], Tmin=(941.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=C1O[CH][C]=CO1(23816)',
    structure = SMILES('C=C1O[CH][C]=CO1'),
    E0 = (210.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07392,0.0264003,4.05872e-05,-6.87865e-08,2.65548e-11,25438.9,15.8225], Tmin=(100,'K'), Tmax=(1017.03,'K')), NASAPolynomial(coeffs=[13.3511,0.0212267,-9.56843e-06,1.96979e-09,-1.49323e-13,21118.8,-48.7356], Tmin=(1017.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'C=C1OC=C1C=O(23817)',
    structure = SMILES('C=C1OC=C1C=O'),
    E0 = (-5.19906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44158,0.0425953,-4.33732e-06,-3.23424e-08,1.69517e-11,-520.937,16.8734], Tmin=(100,'K'), Tmax=(980.314,'K')), NASAPolynomial(coeffs=[17.4067,0.00976463,-3.54435e-06,7.41511e-10,-6.00334e-14,-5203.74,-67.7509], Tmin=(980.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.19906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + ring(Oxetene)"""),
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
    label = '[C]=CO[C]=C(5335)',
    structure = SMILES('[C]=CO[C]=C'),
    E0 = (767.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,357.11,357.921,360.601],'cm^-1')),
        HinderedRotor(inertia=(0.166878,'amu*angstrom^2'), symmetry=1, barrier=(15.3545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.199167,'amu*angstrom^2'), symmetry=1, barrier=(18.1259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22115,0.0373841,-3.89855e-05,2.10503e-08,-4.50258e-12,92393.7,20.9316], Tmin=(100,'K'), Tmax=(1138.69,'K')), NASAPolynomial(coeffs=[9.70072,0.0111097,-4.37411e-06,7.86347e-10,-5.36191e-14,90690.3,-16.1245], Tmin=(1138.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(767.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(C=CJO)"""),
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
    E0 = (282.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (366.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (541.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (479.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (282.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (558.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (395.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (441.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (696.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (375.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (424.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (430.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (734.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (633.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (626.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (470.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (409.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (408.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (339.286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (322.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (579.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (470.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (653.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (911.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (515.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (522.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (445.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (525.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (341.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (291.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (835.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=C=O(598)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=C1OC1[C]=C[O](23793)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(4.8958e+11,'s^-1'), n=-0.055489, Ea=(83.6851,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=[C]OC=C=C=O(23794)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C#COC=C=C[O](23795)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C=O(598)', '[CH]=C=C[O](8556)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.000343471,'m^3/(mol*s)'), n=2.65, Ea=(72.3301,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R_R;CdsJ=Cdd] for rate rule [Od_Cdd;CdsJ=Cdd]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 67.4 to 72.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]OC=C=[C]O(23796)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=COC=C=C[O](23797)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=CO[C]=C=C[O](23798)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.85004e+08,'s^-1'), n=1.20667, Ea=(158.922,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;Cd_H_out_doubleC] + [R3H_SS;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_double;Cd_H_out_doubleC]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]O[C]=C=CO(23799)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=CO[CH][C]=C=O(23800)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.08955e+06,'s^-1'), n=1.58778, Ea=(92.3497,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_Cd;XH_out] + [R5H;Cd_rad_out;XH_out] for rate rule [R5H;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C#CO[CH][C]=CO(23801)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.2360679775
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][C]=O(601)', '[CH]=C=C[O](8556)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', 'C=[C]O[C]=C=C[O](23802)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(8)', 'C=[C]O[CH][C]=C=O(23803)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', 'C#CO[CH][C]=C[O](23804)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=[C]OC=C1[CH]O1(23805)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=C1O[CH]C1=C[O](23806)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=[C]OC1[C]=CO1(23807)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=C1OC=[C]C1[O](23808)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.01659e+11,'s^-1'), n=0.0800815, Ea=(56.4798,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5_linear;doublebond_intra;radadd_intra_cddouble] + [R5_linear;doublebond_intra_CdCdd;radadd_intra] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=COC=C=C=O(23809)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.85329e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]O[C]C=C[O](23810)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]O[CH]C=[C][O](23811)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]=C(584)', '[O]C=C=C[O](22349)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(T)(63)', '[CH]=C=CO[C]=C(19367)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', 'C=[C]OC#CC=O(23812)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]O[C]=CC=O(23813)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=[C]O[CH]C=C=O(23814)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C#CO[CH]C=C[O](23815)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.09763e+09,'s^-1'), n=1.36198, Ea=(242.913,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleDe;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=C1O[CH][C]=CO1(23816)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.19745e+11,'s^-1'), n=0.440371, Ea=(58.418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6;multiplebond_intra;radadd_intra_cddouble] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.2360679775
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]OC=C=C[O](22462)'],
    products = ['C=C1OC=C1C=O(23817)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;CdsinglepriDe_rad_out]
Euclidian distance = 2.2360679775
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=O(373)', '[C]=CO[C]=C(5335)'],
    products = ['C=[C]OC=C=C[O](22462)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4759',
    isomers = [
        'C=[C]OC=C=C[O](22462)',
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
    label = '4759',
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

