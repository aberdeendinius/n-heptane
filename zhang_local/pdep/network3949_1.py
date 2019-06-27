species(
    label = '[CH]=[C]CC[CH][O](13806)',
    structure = SMILES('[CH]=[C]CC[CH][O]'),
    E0 = (690.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,3025,407.5,1350,352.5,220.779,221.091,1615.47],'cm^-1')),
        HinderedRotor(inertia=(0.165671,'amu*angstrom^2'), symmetry=1, barrier=(5.73245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165449,'amu*angstrom^2'), symmetry=1, barrier=(5.73347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165297,'amu*angstrom^2'), symmetry=1, barrier=(5.73192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30193,0.068525,-0.000112272,1.06047e-07,-3.86782e-11,83173.7,25.9957], Tmin=(100,'K'), Tmax=(849.891,'K')), NASAPolynomial(coeffs=[4.26736,0.0344454,-1.66087e-05,3.14835e-09,-2.15093e-13,83396.4,16.4471], Tmin=(849.891,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S)"""),
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
    label = 'C#C[CH2](17441)',
    structure = SMILES('C#C[CH2]'),
    E0 = (328.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2175,525,1131.03,1132.16,1135.9],'cm^-1')),
        HinderedRotor(inertia=(0.154206,'amu*angstrom^2'), symmetry=1, barrier=(3.5455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32026,0.0108736,8.62061e-06,-1.82973e-08,7.68649e-12,39535.3,8.27851], Tmin=(100,'K'), Tmax=(960.555,'K')), NASAPolynomial(coeffs=[6.38511,0.00814486,-2.78734e-06,4.95348e-10,-3.50148e-14,38483.6,-8.79383], Tmin=(960.555,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl)"""),
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
    label = '[CH]=[C]CC=C[O](13802)',
    structure = SMILES('[CH]=[C]CC=C[O]'),
    E0 = (506.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.91324,'amu*angstrom^2'), symmetry=1, barrier=(20.9972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.908941,'amu*angstrom^2'), symmetry=1, barrier=(20.8983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50942,0.0449231,-1.85543e-05,-1.50806e-08,1.09627e-11,60999.4,23.909], Tmin=(100,'K'), Tmax=(952.774,'K')), NASAPolynomial(coeffs=[15.1607,0.011459,-3.41436e-06,5.95962e-10,-4.37422e-14,57315.7,-46.9703], Tmin=(952.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(C=COJ)"""),
)

species(
    label = '[CH][C]=CCC=O(19361)',
    structure = SMILES('[CH][C]=CCC=O'),
    E0 = (479.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,414.087,414.09,414.09,414.091],'cm^-1')),
        HinderedRotor(inertia=(0.407962,'amu*angstrom^2'), symmetry=1, barrier=(49.6403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407962,'amu*angstrom^2'), symmetry=1, barrier=(49.6403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.407962,'amu*angstrom^2'), symmetry=1, barrier=(49.6403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9363,0.0403679,-1.40347e-05,-3.1601e-09,1.89973e-12,57732.2,22.4915], Tmin=(100,'K'), Tmax=(1426.78,'K')), NASAPolynomial(coeffs=[11.3773,0.0265257,-1.27557e-05,2.44436e-09,-1.69006e-13,53753,-30.9148], Tmin=(1426.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=C(18825)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,228.264,228.889,229.07],'cm^-1')),
        HinderedRotor(inertia=(1.35219,'amu*angstrom^2'), symmetry=1, barrier=(50.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27541,0.0127954,9.49515e-06,-1.56026e-08,5.42938e-12,73954,11.3502], Tmin=(100,'K'), Tmax=(1063.31,'K')), NASAPolynomial(coeffs=[4.18965,0.0168435,-6.77763e-06,1.22218e-09,-8.33556e-14,73336.3,4.89309], Tmin=(1063.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=[C]C[CH]C[O](13807)',
    structure = SMILES('[CH]=[C]C[CH]C[O]'),
    E0 = (710.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,3025,407.5,1350,352.5,236.01,1522.11,1522.6],'cm^-1')),
        HinderedRotor(inertia=(0.00552224,'amu*angstrom^2'), symmetry=1, barrier=(9.07354,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227438,'amu*angstrom^2'), symmetry=1, barrier=(9.02301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230163,'amu*angstrom^2'), symmetry=1, barrier=(9.04424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75877,0.0546538,-7.27674e-05,6.48538e-08,-2.40287e-11,85518.8,27.3038], Tmin=(100,'K'), Tmax=(798.712,'K')), NASAPolynomial(coeffs=[4.18654,0.0333982,-1.57642e-05,3.01433e-09,-2.09212e-13,85421.1,17.9534], Tmin=(798.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(710.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CCJCO) + radical(Cds_S) + radical(CCOJ)"""),
)

species(
    label = '[CH]C=CC[CH][O](13805)',
    structure = SMILES('[CH]C=CC[CH][O]'),
    E0 = (564.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,516.977,516.981,516.981,516.982,516.991,516.997],'cm^-1')),
        HinderedRotor(inertia=(0.284387,'amu*angstrom^2'), symmetry=1, barrier=(53.9361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284381,'amu*angstrom^2'), symmetry=1, barrier=(53.9359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284369,'amu*angstrom^2'), symmetry=1, barrier=(53.9359,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53339,0.0596027,-7.16595e-05,6.18407e-08,-2.29724e-11,67954.8,24.9606], Tmin=(100,'K'), Tmax=(793.18,'K')), NASAPolynomial(coeffs=[3.63161,0.0403259,-1.87606e-05,3.55792e-09,-2.46098e-13,67895.5,17.0484], Tmin=(793.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=CCC[O](19537)',
    structure = SMILES('[CH][C]=CCC[O]'),
    E0 = (621.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.20005,0.0403643,-1.8436e-05,3.33616e-09,-2.10758e-13,74796.3,19.7544], Tmin=(100,'K'), Tmax=(2951.51,'K')), NASAPolynomial(coeffs=[39.6293,-0.00217202,-2.91521e-07,2.22941e-11,3.48693e-15,50315.4,-200.466], Tmin=(2951.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(CCOJ)"""),
)

species(
    label = '[CH]=CC[CH][CH][O](19538)',
    structure = SMILES('[CH]=CC[CH][CH][O]'),
    E0 = (652.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24282,0.0471354,-1.6265e-05,-6.92304e-08,7.97199e-11,78576.7,24.8082], Tmin=(100,'K'), Tmax=(458.862,'K')), NASAPolynomial(coeffs=[5.91761,0.0308304,-1.43818e-05,2.73635e-09,-1.89433e-13,78073.9,8.1378], Tmin=(458.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCsJOH) + radical(CCOJ) + radical(Cds_P) + radical(CCJCO)"""),
)

species(
    label = '[CH]=[C]C[CH][CH]O(19304)',
    structure = SMILES('[CH]=[C]C[CH][CH]O'),
    E0 = (665.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,3120,650,792.5,1650,3000,3050,390,425,1340,1360,335,370,312.505,312.517],'cm^-1')),
        HinderedRotor(inertia=(0.00172619,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00172616,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104217,'amu*angstrom^2'), symmetry=1, barrier=(7.22254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104213,'amu*angstrom^2'), symmetry=1, barrier=(7.2226,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33098,0.063328,-7.94697e-05,5.16512e-08,-1.1633e-11,80073.9,27.5517], Tmin=(100,'K'), Tmax=(669.316,'K')), NASAPolynomial(coeffs=[9.3046,0.0241411,-1.062e-05,1.97083e-09,-1.34844e-13,78816.9,-9.1316], Tmin=(669.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(CCJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][C]=CC[CH][O](17845)',
    structure = SMILES('[CH2][C]=CC[CH][O]'),
    E0 = (582.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,3025,407.5,1350,352.5,224.676,1713.33,1713.33],'cm^-1')),
        HinderedRotor(inertia=(0.269387,'amu*angstrom^2'), symmetry=1, barrier=(9.64959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269407,'amu*angstrom^2'), symmetry=1, barrier=(9.64961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34394,'amu*angstrom^2'), symmetry=1, barrier=(48.1353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45678,0.0623132,-9.01891e-05,8.0812e-08,-2.92088e-11,70200.4,24.8796], Tmin=(100,'K'), Tmax=(819.946,'K')), NASAPolynomial(coeffs=[5.0229,0.0332781,-1.57819e-05,3.00369e-09,-2.07155e-13,70006.8,10.7688], Tmin=(819.946,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(Cds_S) + radical(CCsJOH)"""),
)

species(
    label = '[CH][C]=CC[CH]O(19305)',
    structure = SMILES('[CH][C]=CC[CH]O'),
    E0 = (576.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2750,2850,1437.5,1250,1305,750,350,3615,1277.5,1000,1685,370,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14669,0.0668163,-8.34861e-05,6.47362e-08,-2.1018e-11,69429.7,25.9636], Tmin=(100,'K'), Tmax=(805.194,'K')), NASAPolynomial(coeffs=[7.24787,0.0332093,-1.4736e-05,2.72743e-09,-1.85933e-13,68554.1,-1.48528], Tmin=(805.194,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(AllylJ2_triplet) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C[CH][CH][O](17844)',
    structure = SMILES('C=[C]C[CH][CH][O]'),
    E0 = (643.617,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,3000,3050,390,425,1340,1360,335,370,194.398,195.907,198.028,1521.32],'cm^-1')),
        HinderedRotor(inertia=(0.00432632,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00448341,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00460007,'amu*angstrom^2'), symmetry=1, barrier=(7.55266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58091,0.0596599,-8.65415e-05,7.86384e-08,-2.87308e-11,77490.1,27.0881], Tmin=(100,'K'), Tmax=(823.391,'K')), NASAPolynomial(coeffs=[4.47847,0.033179,-1.57025e-05,2.98622e-09,-2.05803e-13,77433.5,16.2256], Tmin=(823.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = '[CH]=[C]C[CH][CH][O](19539)',
    structure = SMILES('[CH]=[C]C[CH][CH][O]'),
    E0 = (890.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3120,650,792.5,1650,3000,3050,390,425,1340,1360,335,370,274.396,1447.08,1447.45],'cm^-1')),
        HinderedRotor(inertia=(0.106505,'amu*angstrom^2'), symmetry=1, barrier=(5.70325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106589,'amu*angstrom^2'), symmetry=1, barrier=(5.70084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106514,'amu*angstrom^2'), symmetry=1, barrier=(5.69514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48262,0.06347,-0.000104167,9.71491e-08,-3.5031e-11,107211,27.8438], Tmin=(100,'K'), Tmax=(850.178,'K')), NASAPolynomial(coeffs=[4.88215,0.0300626,-1.45025e-05,2.74782e-09,-1.87621e-13,107262,15.6959], Tmin=(850.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(890.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(Cds_P) + radical(CCOJ) + radical(CCsJOH) + radical(CCJCO)"""),
)

species(
    label = '[CH][C]=CC[CH][O](19540)',
    structure = SMILES('[CH][C]=CC[CH][O]'),
    E0 = (802.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,531.118,531.119,531.119,531.119,531.119,531.12],'cm^-1')),
        HinderedRotor(inertia=(0.271951,'amu*angstrom^2'), symmetry=1, barrier=(54.4379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271951,'amu*angstrom^2'), symmetry=1, barrier=(54.4379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271951,'amu*angstrom^2'), symmetry=1, barrier=(54.4379,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45102,0.0648094,-9.86543e-05,9.40038e-08,-3.51313e-11,96560.2,25.728], Tmin=(100,'K'), Tmax=(844.631,'K')), NASAPolynomial(coeffs=[2.68846,0.0393913,-1.87808e-05,3.54494e-09,-2.42209e-13,97048.8,24.0968], Tmin=(844.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(802.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCsJOH) + radical(CCOJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=CCC=C[O](12776)',
    structure = SMILES('[CH]=CCC=C[O]'),
    E0 = (268.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.964194,'amu*angstrom^2'), symmetry=1, barrier=(22.1687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96535,'amu*angstrom^2'), symmetry=1, barrier=(22.1953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3714.19,'J/mol'), sigma=(6.10713,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=580.15 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54104,0.0403358,6.2145e-06,-4.43749e-08,2.19672e-11,32396.1,23.322], Tmin=(100,'K'), Tmax=(945.635,'K')), NASAPolynomial(coeffs=[16.3545,0.0119522,-3.13312e-06,5.46156e-10,-4.2351e-14,28061.9,-55.4194], Tmin=(945.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C#C[CH]CC[O](19541)',
    structure = SMILES('C#C[CH]CC[O]'),
    E0 = (337.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63858,0.0561205,-7.21531e-05,5.87321e-08,-1.92566e-11,40719.7,22.0443], Tmin=(100,'K'), Tmax=(895.404,'K')), NASAPolynomial(coeffs=[5.42687,0.0297023,-1.19906e-05,2.09572e-09,-1.36966e-13,40421.9,6.312], Tmin=(895.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.889,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(CCOJ)"""),
)

species(
    label = '[CH]=[C]CC([CH2])[O](17830)',
    structure = SMILES('[CH]=[C]CC([CH2])[O]'),
    E0 = (714.115,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3120,650,792.5,1650,3000,3100,440,815,1455,1000,358.676,359.105],'cm^-1')),
        HinderedRotor(inertia=(0.00130755,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134185,'amu*angstrom^2'), symmetry=1, barrier=(12.2634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134345,'amu*angstrom^2'), symmetry=1, barrier=(12.2604,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.28011,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34582,0.0599191,-6.53843e-05,3.89316e-08,-9.4385e-12,85982.6,25.5327], Tmin=(100,'K'), Tmax=(994.707,'K')), NASAPolynomial(coeffs=[10.4569,0.0232814,-1.01362e-05,1.90406e-09,-1.325e-13,84170,-18.3748], Tmin=(994.707,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(714.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Cds_P) + radical(CJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=C)C[CH][O](17853)',
    structure = SMILES('[CH]C(=C)C[CH][O]'),
    E0 = (561.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3025,407.5,1350,352.5,538.001,538.002,538.002,538.002,538.004,538.005],'cm^-1')),
        HinderedRotor(inertia=(0.261591,'amu*angstrom^2'), symmetry=1, barrier=(53.7307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261591,'amu*angstrom^2'), symmetry=1, barrier=(53.7306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261595,'amu*angstrom^2'), symmetry=1, barrier=(53.7306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35292,0.0638255,-8.0873e-05,6.93227e-08,-2.50298e-11,67596.7,24.5592], Tmin=(100,'K'), Tmax=(805.735,'K')), NASAPolynomial(coeffs=[4.47186,0.0392433,-1.81715e-05,3.42872e-09,-2.36097e-13,67389.4,12.0186], Tmin=(805.735,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(561.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C1CCC1[O](19542)',
    structure = SMILES('[CH]=C1CCC1[O]'),
    E0 = (405.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[11.3792,-0.0282107,0.000142573,-1.3138e-07,3.00787e-11,48327.7,-20.5721], Tmin=(100,'K'), Tmax=(1707.26,'K')), NASAPolynomial(coeffs=[71.3562,0.0224526,-6.99154e-05,1.71875e-08,-1.28148e-12,-14.2696,-423.608], Tmin=(1707.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(405.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Cds_P) + radical(CC(C)OJ)"""),
)

species(
    label = '[CH]=[C]CCC=O(13813)',
    structure = SMILES('[CH]=[C]CCC=O'),
    E0 = (361.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3120,650,792.5,1650,2782.5,750,1395,475,1775,1000,254.36],'cm^-1')),
        HinderedRotor(inertia=(0.155474,'amu*angstrom^2'), symmetry=1, barrier=(7.13828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155476,'amu*angstrom^2'), symmetry=1, barrier=(7.13828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155472,'amu*angstrom^2'), symmetry=1, barrier=(7.13828,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24495,0.0455976,-2.28028e-05,-3.4221e-08,4.34398e-11,43575,21.6415], Tmin=(100,'K'), Tmax=(487.529,'K')), NASAPolynomial(coeffs=[5.34579,0.0315735,-1.4782e-05,2.84605e-09,-1.99662e-13,43137,7.51782], Tmin=(487.529,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH][O](751)',
    structure = SMILES('[CH][O]'),
    E0 = (424.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([815.726,815.726,3402.81],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86392,-0.000399472,1.49306e-05,-2.12194e-08,8.78636e-12,51105.5,7.21901], Tmin=(100,'K'), Tmax=(905.857,'K')), NASAPolynomial(coeffs=[5.97079,-0.000856178,1.03779e-06,-2.14004e-10,1.3909e-14,50360.8,-4.74054], Tmin=(905.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-OsHHH) + radical(H3COJ) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = '[CH]=[C]C[CH2](5188)',
    structure = SMILES('[CH]=[C]C[CH2]'),
    E0 = (673.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.0131566,'amu*angstrom^2'), symmetry=1, barrier=(6.38658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429434,'amu*angstrom^2'), symmetry=1, barrier=(9.87353,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.6422,0.0283978,-1.83267e-05,6.21843e-09,-8.86535e-13,81072.8,17.6148], Tmin=(100,'K'), Tmax=(1566.78,'K')), NASAPolynomial(coeffs=[7.52451,0.015933,-6.39293e-06,1.14053e-09,-7.62755e-14,79542.9,-8.13177], Tmin=(1566.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(673.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C[CH][O](430)',
    structure = SMILES('[CH2]C[CH][O]'),
    E0 = (330.898,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1568.04,1573.36],'cm^-1')),
        HinderedRotor(inertia=(0.245918,'amu*angstrom^2'), symmetry=1, barrier=(5.65413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.245609,'amu*angstrom^2'), symmetry=1, barrier=(5.64703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54439,0.0379772,-5.80752e-05,5.79428e-08,-2.24013e-11,39844.5,17.4831], Tmin=(100,'K'), Tmax=(844.872,'K')), NASAPolynomial(coeffs=[2.07504,0.0262311,-1.24218e-05,2.35033e-09,-1.60939e-13,40422.3,22.6186], Tmin=(844.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.898,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[C]=[CH](18830)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([562.459,1392.74,3112.83],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86574,0.0023643,1.31489e-06,-2.64796e-09,1.00932e-12,101918,6.62295], Tmin=(100,'K'), Tmax=(1016.96,'K')), NASAPolynomial(coeffs=[4.17915,0.00251347,-9.43473e-07,1.6872e-10,-1.15831e-14,101782,4.75427], Tmin=(1016.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC[C][O](19543)',
    structure = SMILES('[CH]=[C]CC[C][O]'),
    E0 = (971.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,283.924,283.924,283.924,4000],'cm^-1')),
        HinderedRotor(inertia=(0.545418,'amu*angstrom^2'), symmetry=1, barrier=(31.2002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159461,'amu*angstrom^2'), symmetry=1, barrier=(9.1219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159461,'amu*angstrom^2'), symmetry=1, barrier=(9.12192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34343,0.0659121,-0.000106872,9.79121e-08,-3.51541e-11,116933,25.0294], Tmin=(100,'K'), Tmax=(829.345,'K')), NASAPolynomial(coeffs=[5.98434,0.0291814,-1.44901e-05,2.79185e-09,-1.92974e-13,116657,6.48278], Tmin=(829.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(971.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CH2_triplet) + radical(Cds_P) + radical(CCOJ)"""),
)

species(
    label = '[C]=[C]CC[CH][O](19544)',
    structure = SMILES('[C]=[C]CC[CH][O]'),
    E0 = (1001.82,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,180,180,180,1979.42],'cm^-1')),
        HinderedRotor(inertia=(0.158348,'amu*angstrom^2'), symmetry=1, barrier=(3.64074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159641,'amu*angstrom^2'), symmetry=1, barrier=(3.67046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158998,'amu*angstrom^2'), symmetry=1, barrier=(3.65567,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32434,0.0715619,-0.00013349,1.33041e-07,-4.93763e-11,120575,25.8894], Tmin=(100,'K'), Tmax=(864.459,'K')), NASAPolynomial(coeffs=[2.77554,0.0347749,-1.74772e-05,3.33081e-09,-2.26585e-13,121448,25.5986], Tmin=(864.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1001.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCsJOH) + radical(CCOJ) + radical(CdCdJ2_triplet)"""),
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
    E0 = (690.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (725.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (702.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (690.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (716.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (847.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (881.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (806.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (832.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (846.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (835.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (814.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (899.826,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (975.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1102.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1013.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (754.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (754.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (871.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (936.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (699.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (690.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1132.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1212.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1183.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1213.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['C=C[O](594)', 'C#C[CH2](17441)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['H(8)', '[CH]=[C]CC=C[O](13802)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', '[CH][C]=CCC=O(19361)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C[O](594)', '[CH][C]=C(18825)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(0.123469,'m^3/(mol*s)'), n=2.00579, Ea=(101.342,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 98.5 to 101.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH][O](719)', 'C#C[CH2](17441)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(0.00755369,'m^3/(mol*s)'), n=2.43138, Ea=(27.0165,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]C[CH]C[O](13807)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH]C=CC[CH][O](13805)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(9.57388e+10,'s^-1'), n=0.905639, Ea=(190.563,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH][C]=CCC[O](19537)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH]=CC[CH][CH][O](19538)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH]=[C]C[CH][CH]O(19304)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.30814e+09,'s^-1'), n=1.19923, Ea=(155.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH2][C]=CC[CH][O](17845)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH][C]=CC[CH]O(19305)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.04154e+06,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['C=[C]C[CH][CH][O](17844)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][O](719)', '[CH][C]=C(18825)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.9789e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', '[CH]=[C]C[CH][CH][O](19539)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(8)', '[CH][C]=CC[CH][O](19540)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH]=CCC=C[O](12776)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['C#C[CH]CC[O](19541)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]CC([CH2])[O](17830)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH]C(=C)C[CH][O](17853)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH]=C1CCC1[O](19542)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=[C]CC[CH][O](13806)'],
    products = ['[CH]=[C]CCC=O(13813)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH][O](751)', '[CH]=[C]C[CH2](5188)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C[CH][O](430)', '[C]=[CH](18830)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(8)', '[CH]=[C]CC[C][O](19543)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(8)', '[C]=[C]CC[CH][O](19544)'],
    products = ['[CH]=[C]CC[CH][O](13806)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

network(
    label = '3949',
    isomers = [
        '[CH]=[C]CC[CH][O](13806)',
    ],
    reactants = [
        ('C=C[O](594)', 'C#C[CH2](17441)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '3949',
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

