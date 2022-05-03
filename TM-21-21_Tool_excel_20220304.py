# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:49:18 2021

@author: Jianchuan Tan at PNNL (tanj239)
"""
#--------------------------- Import Packages ---------------------------------#
import numpy as np
import pandas as pd
import math
import sys

#------------------------------- Excel Input ---------------------------------#
df = pd.read_excel('Data_E21.xlsx')
NumTest = len([item for item in list(df.loc[0]) if math.isnan(item)==False])

#-------------------------------- JSON Input ---------------------------------#
fluxData = {
    "condition1": {
        "testData": {
            "caseTemperature": df[df.columns[0]][0],
            "driveCurrent": df[df.columns[0]][1],
            "numberOfUnitsTested": df[df.columns[0]][2],
            "numberOfFailures": df[df.columns[0]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[0]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[1]][5::]) if math.isnan(item)==False]
        }
    },
     "condition2": {
        "testData": {
            "caseTemperature": df[df.columns[2]][0],
            "driveCurrent": df[df.columns[2]][1],
            "numberOfUnitsTested": df[df.columns[2]][2],
            "numberOfFailures": df[df.columns[2]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[2]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[3]][5::]) if math.isnan(item)==False]
        }
    },
     "condition3": {
        "testData": {
            "caseTemperature": df[df.columns[4]][0],
            "driveCurrent": df[df.columns[4]][1],
            "numberOfUnitsTested": df[df.columns[4]][2],
            "numberOfFailures": df[df.columns[4]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[4]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[5]][5::]) if math.isnan(item)==False]
        }
    },
     "condition4": {
        "testData": {
            "caseTemperature": df[df.columns[6]][0],
            "driveCurrent": df[df.columns[6]][1],
            "numberOfUnitsTested": df[df.columns[6]][2],
            "numberOfFailures": df[df.columns[6]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[6]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[7]][5::]) if math.isnan(item)==False]
        }
    },
     "condition5": {
        "testData": {
            "caseTemperature": df[df.columns[8]][0],
            "driveCurrent": df[df.columns[8]][1],
            "numberOfUnitsTested": df[df.columns[8]][2],
            "numberOfFailures": df[df.columns[8]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[8]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[9]][5::]) if math.isnan(item)==False]
        }
    },
     "condition6": {
        "testData": {
            "caseTemperature": df[df.columns[10]][0],
            "driveCurrent": df[df.columns[10]][1],
            "numberOfUnitsTested": df[df.columns[10]][2],
            "numberOfFailures": df[df.columns[10]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[10]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[11]][5::]) if math.isnan(item)==False]
        }
    },
     "condition7": {
        "testData": {
            "caseTemperature": df[df.columns[12]][0],
            "driveCurrent": df[df.columns[12]][1],
            "numberOfUnitsTested": df[df.columns[12]][2],
            "numberOfFailures": df[df.columns[12]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[12]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[13]][5::]) if math.isnan(item)==False]
        }
    },
     "condition8": {
        "testData": {
            "caseTemperature": df[df.columns[14]][0],
            "driveCurrent": df[df.columns[14]][1],
            "numberOfUnitsTested": df[df.columns[14]][2],
            "numberOfFailures": df[df.columns[14]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[14]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[15]][5::]) if math.isnan(item)==False]
        }
    },
     "condition9": {
        "testData": {
            "caseTemperature": df[df.columns[16]][0],
            "driveCurrent": df[df.columns[16]][1],
            "numberOfUnitsTested": df[df.columns[16]][2],
            "numberOfFailures": df[df.columns[16]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[16]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[17]][5::]) if math.isnan(item)==False]
        }
    },
     "condition10": {
        "testData": {
            "caseTemperature": df[df.columns[18]][0],
            "driveCurrent": df[df.columns[18]][1],
            "numberOfUnitsTested": df[df.columns[18]][2],
            "numberOfFailures": df[df.columns[18]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[18]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[19]][5::]) if math.isnan(item)==False]
        }
    },
     "condition11": {
        "testData": {
            "caseTemperature": df[df.columns[20]][0],
            "driveCurrent": df[df.columns[20]][1],
            "numberOfUnitsTested": df[df.columns[20]][2],
            "numberOfFailures": df[df.columns[20]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[20]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[21]][5::]) if math.isnan(item)==False]
        }
    },
     "condition12": {
        "testData": {
            "caseTemperature": df[df.columns[22]][0],
            "driveCurrent": df[df.columns[22]][1],
            "numberOfUnitsTested": df[df.columns[22]][2],
            "numberOfFailures": df[df.columns[22]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[22]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[23]][5::]) if math.isnan(item)==False]
        }
    },
     "condition13": {
        "testData": {
            "caseTemperature": df[df.columns[24]][0],
            "driveCurrent": df[df.columns[24]][1],
            "numberOfUnitsTested": df[df.columns[24]][2],
            "numberOfFailures": df[df.columns[24]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[24]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[25]][5::]) if math.isnan(item)==False]
        }
    },
     "condition14": {
        "testData": {
            "caseTemperature": df[df.columns[26]][0],
            "driveCurrent": df[df.columns[26]][1],
            "numberOfUnitsTested": df[df.columns[26]][2],
            "numberOfFailures": df[df.columns[26]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[26]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[27]][5::]) if math.isnan(item)==False]
        }
    },
     "condition15": {
        "testData": {
            "caseTemperature": df[df.columns[28]][0],
            "driveCurrent": df[df.columns[28]][1],
            "numberOfUnitsTested": df[df.columns[28]][2],
            "numberOfFailures": df[df.columns[28]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[28]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[29]][5::]) if math.isnan(item)==False]
        }
    },
     "condition16": {
        "testData": {
            "caseTemperature": df[df.columns[30]][0],
            "driveCurrent": df[df.columns[30]][1],
            "numberOfUnitsTested": df[df.columns[30]][2],
            "numberOfFailures": df[df.columns[30]][3]
        },
        "measuredData": {
            "hours": [item for item in list(df[df.columns[30]][5::]) if math.isnan(item)==False],
            "flux": [item for item in list(df[df.columns[31]][5::]) if math.isnan(item)==False]
        }
    }
}

interpolationData = {
    "interpolationParameters": {
        "inSituTemperature": 117,
        "inSituCurrent": 700,
        "fluxType": "Luminous"
    },
    "ledPackageCharacterstics":{
        "ratedMaxCurrent": 1000,
        "nominalForwardVoltage": 3.2
    }
}

reportedValues = {
    "fluxMaintenanceLevel(%)": 80,
    "projectedFluxMaintenanceHours": 35000
}

#------------------------------ Build DataFrame ------------------------------#
CaseTempC_raw = []
CaseTempC = []
Current_raw = []
Current = []
NumUnits = []
NumFail = []
Hours = []
Flux = []
Hours_ori = []
Flux_ori = []
BoldRed = []
for i in [x for x in range(1, NumTest+1)]:
    CaseTempC_raw.append(fluxData["condition"+str(i)]["testData"]["caseTemperature"])
    CaseTempC.append(fluxData["condition"+str(i)]["testData"]["caseTemperature"])
    Current_raw.append(fluxData["condition"+str(i)]["testData"]["driveCurrent"])
    Current.append(fluxData["condition"+str(i)]["testData"]["driveCurrent"])
    NumUnits.append(fluxData["condition"+str(i)]["testData"]["numberOfUnitsTested"])
    NumFail.append(fluxData["condition"+str(i)]["testData"]["numberOfFailures"])
    Hours.append(fluxData["condition"+str(i)]["measuredData"]["hours"])
    Flux.append(fluxData["condition"+str(i)]["measuredData"]["flux"])
    Hours_ori.append(fluxData["condition"+str(i)]["measuredData"]["hours"])
    Flux_ori.append(fluxData["condition"+str(i)]["measuredData"]["flux"])

maxLEDpower = interpolationData["ledPackageCharacterstics"]["nominalForwardVoltage"] * \
            interpolationData["ledPackageCharacterstics"]["ratedMaxCurrent"] / 1000
            
if interpolationData["interpolationParameters"]["fluxType"] == "Luminous":
    Type = "L"
elif interpolationData["interpolationParameters"]["fluxType"] == "Photon":
    Type = "Q"
elif interpolationData["interpolationParameters"]["fluxType"] == "Radiant":
    Type = "R"
else:
    sys.exit("Error: Flux Type must be Luminous, or Photon, or Radiant.")
    
#---------------- Group, Average, and Replace Case Temperatures --------------#
TempUncertainty = 1.5
TempGroupIndex = []

TempGroupIndex.append(1)
for i in range(1, len(CaseTempC_raw)):
    TempGroupIndex.append(0)
    for j in range(i):
        if abs(CaseTempC_raw[i] - CaseTempC_raw[j]) < TempUncertainty:
            TempGroupIndex[i] = TempGroupIndex[j]
        else:
            continue
    if TempGroupIndex[i] == 0:
        TempGroupIndex[i] = max(TempGroupIndex) + 1
    else:
        continue

for idx_temp in range(1, max(TempGroupIndex)+1):
    idx = [i for i,x in enumerate(TempGroupIndex) if x==idx_temp]
    sum_temp = 0
    mean_temp = 0
    for j in idx:
        sum_temp = sum_temp + CaseTempC_raw[j]
    mean_temp = round(sum_temp / len(idx), 1)
    for j in idx:
        CaseTempC[j] = mean_temp
        
CaseTempK = [x+273.15 for x in CaseTempC]

#-------------------- Group, Average, and Replace Current --------------------#
CurrentUncertainty = 5
CurrentGroupIndex = []

CurrentGroupIndex.append(1)
for i in range(1, len(Current_raw)):
    CurrentGroupIndex.append(0)
    for j in range(i):
        if abs(Current_raw[i] - Current_raw[j]) < CurrentUncertainty:
            CurrentGroupIndex[i] = CurrentGroupIndex[j]
        else:
            continue
    if CurrentGroupIndex[i] == 0:
        CurrentGroupIndex[i] = max(CurrentGroupIndex) + 1
    else:
        continue

for idx_curr in range(1, max(CurrentGroupIndex)+1):
    idx = [i for i,x in enumerate(CurrentGroupIndex) if x==idx_curr]
    sum_curr = 0
    mean_curr = 0
    for j in idx:
        sum_curr = sum_curr + Current_raw[j]
    mean_curr = round(sum_curr / len(idx))
    for j in idx:
        Current[j] = mean_curr

#--------------------------- Apply Interpolation Rules -----------------------#
TestTempC = interpolationData["interpolationParameters"]["inSituTemperature"]
TestCurrent = interpolationData["interpolationParameters"]["inSituCurrent"]

if TestTempC < min(CaseTempC):
    TestTempC = min(CaseTempC)

if TestCurrent < min(Current):
    TestCurrent = min(Current)

# elif TestCurrent < 350:
#     TestCurrent = 350
#     break
# elif TestCurrent > 1000 and TestTempC <= 85:
#     print("Error: Test Current should not exceed 1000 mA at this Case Temperature. \nPlease enter again.")
# elif TestCurrent > 1200 and 85 < TestTempC <= 105:
#     print("Error: Test Current should not exceed 1200 mA at this Case Temperature. \nPlease enter again.")
# elif TestCurrent > 700 and 105 < TestTempC <= 120:
#     print("Error: Test Current should not exceed 700 mA at this Case Temperature. \nPlease enter again.")
# else:
#     break

#---------------------------- Projection Rules -------------------------------#
cond_range = range(len(CaseTempC))
proj_limit = []
for i in cond_range:
    #Hours[i].dropna(inplace=True)
    #Flux[i].dropna(inplace=True)
    if max(Hours[i]) < 5952:
        print("Error: Condition "+ str(i+1) +" test duration must be >= 5952 hrs.")
        break
        sys.exit()
    
    if (NumUnits[i] - NumFail[i]) >= 20:
        limit = 6*max(Hours[i])
    elif (NumUnits[i] - NumFail[i]) <=19 and (NumUnits[i] - NumFail[i])>=10:
        limit = 5.5*max(Hours[i])
    else:
        print("Error: Condition "+ str(i+1) +" DUT quantity must be >= 10.")
        break
        sys.exit()
    proj_limit.append(limit)

#------------------ Drive Current Interpolation Restrictions -----------------#
if maxLEDpower > 0.6:
    Current_Ratio_Limit = 2.0
elif maxLEDpower <= 0.6:
    Current_Ratio_Limit = 1.5

#---------------------------- Data Selection for Modeling --------------------#
for i in cond_range:
    cut_idx = 0
    if max(Hours[i]) >= 5952 and max(Hours[i]) <= 10048:
        cut_idx = next(x for x, val in enumerate(Hours[i]) if val > max(Hours[i])-5000) - 1
        for j in range(cut_idx):
            Hours[i].pop(0)
            Flux[i].pop(0)
    elif max(Hours[i]) > 10048:
        cut_idx = next(x for x, val in enumerate(Hours[i]) if val > max(Hours[i])/2) - 1
        for j in range(cut_idx):
            Hours[i].pop(0)
            Flux[i].pop(0)
    for j in range(len(Hours[i])-1):
        if (abs(Hours[i][j+1] - Hours[i][j]) > 1048):
            sys.exit("Error: Condition "+ str(i+1) +" has interval larger than 1048 hrs.")
            break

#---------------------------- Linear Fit Calculation -------------------------#
Flux_log = []
m = []
b = []
B = []
alpha = []
Lxx_calc = []
Lxx_report = []
for i in cond_range:
    Fluxi_log = []
    for j in range(len(Hours[i])):
        Fluxi_log.append(np.log(Flux[i][j]))
    Flux_log.append(Fluxi_log)
    
    xi = list(Hours[i])
    yi = Flux_log[i]
    xyi = [xj * yj for xj, yj in zip(xi, yi)]
    xxi = [xj * xj for xj in xi]
    ni = len(Hours[i])
    mi = (ni*sum(xyi) - sum(xi)*sum(yi)) / (ni*sum(xxi) - sum(xi)**2)
    bi = (sum(yi) - mi*sum(xi)) / ni
    m.append(mi)
    b.append(bi)
    
    Bi = np.exp(bi)
    alphai = -1*mi
    if alphai < 2e-6:   # Apply minimum alpha rule
        alphai_new = 2e-6
        Bi_new = Flux[i][len(Flux[i])-1] * np.exp(2e-6 * max(Hours[i]))
        boldred = 1
    else:
        alphai_new = alphai
        Bi_new = Bi
        boldred = 0
    B.append(Bi_new)
    alpha.append(alphai_new)
    BoldRed.append(boldred)
    
    Lxxi_calc = np.log(Bi_new / (reportedValues["fluxMaintenanceLevel(%)"] / 100)) / alphai_new
    if Lxxi_calc <= proj_limit[i]:
        Lxxi_report = "= " + str("{:,}".format(int(round(Lxxi_calc/100))*100))
    elif Lxxi_calc > proj_limit[i]:
        Lxxi_report = "> " + str("{:,}".format(int(proj_limit[i])))
    Lxx_calc.append(Lxxi_calc)
    Lxx_report.append(Lxxi_report)

####################### Interpolation Between Test Conditions #################
#------------------------ Find conditions for interpolation ------------------#
CaseTempC_sort_set = sorted(set(CaseTempC))
Current_sort_set = sorted(set(Current))

temp_high = [x for x in CaseTempC_sort_set if x >= TestTempC][0]
temp_low = [x for x in CaseTempC_sort_set if x <= TestTempC][-1]
curr_high = [x for x in Current_sort_set if x >= TestCurrent][0]
curr_low = [x for x in Current_sort_set if x <= TestCurrent][-1]

idx_temp_high = [x for x, val in enumerate(CaseTempC) if val==temp_high]
idx_temp_low = [x for x, val in enumerate(CaseTempC) if val==temp_low]
idx_curr_high = [x for x, val in enumerate(Current) if val==curr_high]
idx_curr_low = [x for x, val in enumerate(Current) if val==curr_low]

idx_tempH_currH = list(set(idx_temp_high).intersection(idx_curr_high))[0]
idx_tempH_currL = list(set(idx_temp_high).intersection(idx_curr_low))[0]
idx_tempL_currH = list(set(idx_temp_low).intersection(idx_curr_high))[0]
idx_tempL_currL = list(set(idx_temp_low).intersection(idx_curr_low))[0]

ProjLimit = min([proj_limit[q] for q in [idx_tempL_currL, idx_tempL_currH, idx_tempH_currL, idx_tempH_currH]])

if TestTempC > CaseTempC[idx_tempH_currH] or TestCurrent > Current[idx_tempH_currH]:
    print("Error: Test condition is out of Interpolation Rigion." \
          "\nPlease enter again.")
    sys.exit()

#---------------------------- Current Interpolation --------------------------#
# At lower temperature
if idx_tempL_currH == idx_tempL_currL:
    m_alpha_tempL = 0
    b_alpha_tempL = alpha[idx_tempL_currL]
    alpha_interp_tempL = b_alpha_tempL
    m_B_tempL = 0
    b_B_tempL = B[idx_tempL_currL]
    B_interp_tempL = b_B_tempL
else:    
    m_alpha_tempL = (alpha[idx_tempL_currH] - alpha[idx_tempL_currL]) / (Current[idx_tempL_currH] - Current[idx_tempL_currL])
    b_alpha_tempL = alpha[idx_tempL_currL] - m_alpha_tempL * Current[idx_tempL_currL]
    alpha_interp_tempL = m_alpha_tempL * TestCurrent + b_alpha_tempL
    m_B_tempL = (B[idx_tempL_currH] - B[idx_tempL_currL]) / (Current[idx_tempL_currH] - Current[idx_tempL_currL])
    b_B_tempL = B[idx_tempL_currL] - m_B_tempL * Current[idx_tempL_currL]
    B_interp_tempL = m_B_tempL * TestCurrent + b_B_tempL

if Current[idx_tempL_currH]/Current[idx_tempL_currL] > Current_Ratio_Limit or \
    Current[idx_tempH_currH]/Current[idx_tempH_currL] > Current_Ratio_Limit:
        print("Error: Currents for interpolation exceed 'Current Ratio Limit': " + str(Current_Ratio_Limit) + ".")
        sys.exit()

if alpha_interp_tempL == 2e-6:
    boldred_interp_tempL = 1
else:
    boldred_interp_tempL = 0

Lxx_interp_tempL_calc = np.log(B_interp_tempL / (reportedValues["fluxMaintenanceLevel(%)"] / 100)) / alpha_interp_tempL
if Lxx_interp_tempL_calc <= min(proj_limit[idx_tempL_currL], proj_limit[idx_tempL_currH]):
    Lxx_interp_tempL_report = str("{:,}".format(int(round(Lxx_interp_tempL_calc/100))*100))
else:
    Lxx_interp_tempL_report = "> " + str("{:,}".format(int(min(proj_limit[idx_tempL_currL], proj_limit[idx_tempL_currH]))))

# At higher temperature
if idx_tempH_currH == idx_tempH_currL:
    m_alpha_tempH = 0
    b_alpha_tempH = alpha[idx_tempH_currL]
    alpha_interp_tempH = b_alpha_tempH
    m_B_tempH = 0
    b_B_tempH = B[idx_tempH_currL]
    B_interp_tempH = b_B_tempH
else:
    m_alpha_tempH = (alpha[idx_tempH_currH] - alpha[idx_tempH_currL]) / (Current[idx_tempH_currH] - Current[idx_tempH_currL])
    b_alpha_tempH = alpha[idx_tempH_currL] - m_alpha_tempH * Current[idx_tempH_currL]
    alpha_interp_tempH = m_alpha_tempH * TestCurrent + b_alpha_tempH
    m_B_tempH = (B[idx_tempH_currH] - B[idx_tempH_currL]) / (Current[idx_tempH_currH] - Current[idx_tempH_currL])
    b_B_tempH = B[idx_tempH_currL] - m_B_tempH * Current[idx_tempH_currL]
    B_interp_tempH = m_B_tempH * TestCurrent + b_B_tempH

if alpha_interp_tempH == 2e-6:
    boldred_interp_tempH = 1
else:
    boldred_interp_tempH = 0

Lxx_interp_tempH_calc = np.log(B_interp_tempH / (reportedValues["fluxMaintenanceLevel(%)"] / 100)) / alpha_interp_tempH
if Lxx_interp_tempH_calc <= min(proj_limit[idx_tempH_currL], proj_limit[idx_tempH_currH]):
    Lxx_interp_tempH_report = str("{:,}".format(int(round(Lxx_interp_tempH_calc/100))*100))
else:
    Lxx_interp_tempH_report = "> " + str("{:,}".format(int(min(proj_limit[idx_tempH_currL], proj_limit[idx_tempH_currH]))))

#------------------------- Temperature Interpolation -------------------------#
if CaseTempC[idx_tempL_currL] == CaseTempC[idx_tempH_currL]:
    Ea_kB = np.nan
    A = np.nan
    m_B = np.nan
    b_B = np.nan
    alpha_final = alpha_interp_tempL
    B_final = B_interp_tempL
else:
    Ea_kB = (np.log(alpha_interp_tempL) - np.log(alpha_interp_tempH)) / \
        (1/CaseTempK[idx_tempH_currL] - 1/CaseTempK[idx_tempL_currL])
    A = alpha_interp_tempL * np.exp(Ea_kB / CaseTempK[idx_tempL_currL])
    m_B = (B_interp_tempH - B_interp_tempL) / (CaseTempC[idx_tempH_currL] - CaseTempC[idx_tempL_currL])
    b_B = B_interp_tempL - m_B * CaseTempC[idx_tempL_currL]
    alpha_final = A * np.exp(-1*Ea_kB / (TestTempC + 273.15))
    B_final = m_B * TestTempC + b_B

if alpha_final == 2e-6:
    boldred_final = 1
else:
    boldred_final = 0

############################## Calculate Results ##############################
perc = reportedValues["fluxMaintenanceLevel(%)"]
Lxx_final_calc = np.log(B_final / (perc/100)) / alpha_final
if Lxx_final_calc <= ProjLimit:
    Lxx_final_report = str("{:,}".format(int(round(Lxx_final_calc/100))*100))
else:
    Lxx_final_report = "> " + str("{:,}".format(int(ProjLimit)))

ProjFluxMaintenance = B_final * np.exp(-1*alpha_final * reportedValues["projectedFluxMaintenanceHours"])

######################### Output Calculated Results ###########################
print('Test Temperature (degC): ' + str(TestTempC))
print('Test Current (mA): ' + str(TestCurrent))
print('Flux Maintenance Level: ' + str(perc) + '%')
print(Type + str(perc) +" (hours): "+ Lxx_final_report)
print('LED Package Rated Max Current (mA): ' + str(interpolationData["ledPackageCharacterstics"]["ratedMaxCurrent"]))
print('LED package Nominal Vf (V): ' + str(interpolationData["ledPackageCharacterstics"]["nominalForwardVoltage"]))
print('LED package Max Input Power (W): ' + str(maxLEDpower))
print('Projected Flux Maintenance Hours: ' + str("{:,}".format(reportedValues["projectedFluxMaintenanceHours"])))
print('Projected Luminous Flux Maintenance: ' + str("{:.1f}".format(ProjFluxMaintenance * 100)) + '%')
