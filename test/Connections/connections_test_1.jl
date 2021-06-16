#Put your code here Sangchu.  

using CSV, DataFrames

#Import monotonic data from CSV file.
data = CSV.File("/Users/crismoen/.julia/dev/StructuresKit/src/assets/Connections/Zhang_2020_steel_steel_screw_monotonic.csv")

#Interpolation









import XLSX

# Read fitst lab data Excel file. Please change the path if necessary
el = XLSX.readxlsx("D:/JHU Spring 2021/Fastener Testing Data_Zhidong/Trail1/Central_Summary_20200108_PhaseI.xlsx")

# Relate the experiment serial number with its filename
filenames1 = el["Sheet1!G2:G94"]
fileorders1 = el["Sheet1!A2:A94"]

# Read txt. files by order
for i in 1:size(filenames1,1)
    # Creat two temporary recording current data, one for displacement, the other for force
    disp = []
    force = []

    # Concatenate the file name with orders
    file2read = filenames1[i]
    # Please change the path if necessary
    filepath = "D:/JHU Spring 2021/Fastener Testing Data_Zhidong/Fastener Testing Data_Zhidong/"
    filename = filepath * file2read

    # Read the txt. lab data file
    f = open(filename, "r")
    lines = readlines(f, keep = true)

    # Only pick the displacement column and force column
    for line in lines

        spliteddata = split(line)
        choosedisp = spliteddata[[8],:]
        chooseforce = spliteddata[[10],:]

        # Record them
        append!(disp,choosedisp)
        append!(force,chooseforce)

    end

    close(f)

    # load the displacement data and force data in to DataFrame
    using DataFrames

    df = DataFrame(displacement = disp, forcemag = force)

    using CSV

    outputorder = string(i)
    outputtype = ".csv"
    outputfilename = outputorder * outputtype

    CSV.write(outputfilename, df, header = false)

end

##################################################
#     Same process for the second Excel file     #
##################################################

# Please change the path if necessary
el = XLSX.readxlsx("D:/JHU Spring 2021/Fastener Testing Data_Zhidong/Trail1/Central_Summary_20200108_Phase2&3.xlsx")

filenames2 = el["Sheet1!G2:G64"]

fileorders2 = el["Sheet1!A2:A64"]


for i in 1:size(filenames2,1)

    disp = []
    force = []
    file2read = filenames1[i]
    # Please change the path if necessary
    filepath = "D:/JHU Spring 2021/Fastener Testing Data_Zhidong/Fastener Testing Data_Zhidong/"

    filename = filepath * file2read



    f = open(filename, "r")

    lines = readlines(f, keep = true)

    for line in lines

        spliteddata = split(line)
        choosedisp = spliteddata[[8],:]
        chooseforce = spliteddata[[10],:]
        append!(disp,choosedisp)
        append!(force,chooseforce)
    end

    close(f)

    using DataFrames
    df = DataFrame(displacement = disp, forcemag = force)

    for k in 1:size(disp,1)
        disp[k] = parse(Float64,disp[k])
    end

    for j in 1:size(force,1)
        force[j] = parse(Float64,force[j])
    end

    Maxforce = findmax(force)
    println(Maxforce)
    println(disp[Maxforce[2]])
    println(i)

    using CSV

    outputorder = string(i+88)
    outputtype = ".csv"
    outputfilename = outputorder * outputtype
    CSV.write(outputfilename, df, header = false)

end

#############################
#  Plot the Monotonic data  #
#############################

using DataFrames
using CSV
using Plots

# Load the CSV file. Please change the path if necessary
el = DataFrame(CSV.File("D:/JHU Spring 2021/Fastener Testing Data_Zhidong/Trail1/Mon (2).csv"))

# Convert the already read CSV file into array
df_matrix = Matrix(el)

# test_54_8
test_54_8_t1 = convert(Array, df_matrix[2:4,10])
test_54_8_t2 = convert(Array, df_matrix[2:4,11])
test_54_8_Fy1 = convert(Array, df_matrix[2:4,12])
test_54_8_Fy2 = convert(Array, df_matrix[2:4,13])
test_54_8_Fu1 = convert(Array, df_matrix[2:4,14])
test_54_8_Fu2 = convert(Array, df_matrix[2:4,15])
test_54_8_screw = convert(Array, df_matrix[2:4,16])
test_54_8_thickness = convert(Array, df_matrix[1:3,17])

using Dierckx
test_54_8_t1_int = Spline1D(test_54_8_thickness,test_54_8_t1,k=2)
test_54_8_t2_int = Spline1D(test_54_8_thickness,test_54_8_t2,k=2)
test_54_8_Fy1_int = Spline1D(test_54_8_thickness,test_54_8_Fy1,k=2)
test_54_8_Fy2_int = Spline1D(test_54_8_thickness,test_54_8_Fy2,k=2)
test_54_8_Fu1_int = Spline1D(test_54_8_thickness,test_54_8_Fu1,k=2)
test_54_8_Fu2_int = Spline1D(test_54_8_thickness,test_54_8_Fu2,k=2)


# test_54_10
test_54_10_t1 = convert(Array, df_matrix[5:7,10])
test_54_10_t2 = convert(Array, df_matrix[5:7,11])
test_54_10_Fy1 = convert(Array, df_matrix[5:7,12])
test_54_10_Fy2 = convert(Array, df_matrix[5:7,13])
test_54_10_Fu1 = convert(Array, df_matrix[5:7,14])
test_54_10_Fu2 = convert(Array, df_matrix[5:7,15])
test_54_10_screw = convert(Array, df_matrix[5:7,16])
test_54_10_thickness = convert(Array, df_matrix[1:3,17])

using Dierckx
test_54_10_t1_int = Spline1D(test_54_10_thickness,test_54_10_t1,k=2)
test_54_10_t2_int = Spline1D(test_54_10_thickness,test_54_10_t2,k=2)
test_54_10_Fy1_int = Spline1D(test_54_10_thickness,test_54_10_Fy1,k=2)
test_54_10_Fy2_int = Spline1D(test_54_10_thickness,test_54_10_Fy2,k=2)
test_54_10_Fu1_int = Spline1D(test_54_10_thickness,test_54_10_Fu1,k=2)
test_54_10_Fu2_int = Spline1D(test_54_10_thickness,test_54_10_Fu2,k=2)

# test_97_10
test_97_10_t1 = convert(Array, df_matrix[13:15,10])
test_97_10_t2 = convert(Array, df_matrix[13:15,11])
test_97_10_Fy1 = convert(Array, df_matrix[13:15,12])
test_97_10_Fy2 = convert(Array, df_matrix[13:15,13])
test_97_10_Fu1 = convert(Array, df_matrix[13:15,14])
test_97_10_Fu2 = convert(Array, df_matrix[13:15,15])
test_97_10_screw = convert(Array, df_matrix[13:15,16])
test_97_10_thickness = convert(Array, df_matrix[1:3,17])

using Dierckx
test_97_10_t1_int = Spline1D(test_97_10_thickness,test_97_10_t1,k=2)
test_97_10_t2_int = Spline1D(test_97_10_thickness,test_97_10_t2,k=2)
test_97_10_Fy1_int = Spline1D(test_97_10_thickness,test_97_10_Fy1,k=2)
test_97_10_Fy2_int = Spline1D(test_97_10_thickness,test_97_10_Fy2,k=2)
test_97_10_Fu1_int = Spline1D(test_97_10_thickness,test_97_10_Fu1,k=2)
test_97_10_Fu2_int = Spline1D(test_97_10_thickness,test_97_10_Fu2,k=2)

# test_97_xx_30
test_97_xx_30_t1 = convert(Array, df_matrix[15:17,10])
test_97_xx_30_t2 = convert(Array, df_matrix[15:17,11])
test_97_xx_30_Fy1 = convert(Array, df_matrix[15:17,12])
test_97_xx_30_Fy2 = convert(Array, df_matrix[15:17,13])
test_97_xx_30_Fu1 = convert(Array, df_matrix[15:17,14])
test_97_xx_30_Fu2 = convert(Array, df_matrix[15:17,15])
test_97_xx_30_screw = convert(Array, df_matrix[15:17,16])
test_97_xx_30_thickness = convert(Array, df_matrix[4:6,17])

using Dierckx
test_97_xx_30_int = Spline1D(test_97_xx_30_thickness,test_97_xx_30_t1,k=2)
test_97_xx_30_int = Spline1D(test_97_xx_30_thickness,test_97_xx_30_t2,k=2)
test_97_xx_30_int = Spline1D(test_97_xx_30_thickness,test_97_xx_30_Fy1,k=2)
test_97_xx_30_int = Spline1D(test_97_xx_30_thickness,test_97_xx_30_Fy2,k=2)
test_97_xx_30_int = Spline1D(test_97_xx_30_thickness,test_97_xx_30_Fu1,k=2)
test_97_xx_30_int = Spline1D(test_97_xx_30_thickness,test_97_xx_30_Fu2,k=2)

# Visualize the chosen node

 for i in 1:size(df_matrix,1)

    # Load the 4 node's displacement and force
    F = df_matrix[i,6:9]
    D = df_matrix[i,2:5]

    name = df_matrix[i,1]

    df = DataFrame(force = F, displacement = D)
    plot(D,F)
    Plots.display(plot(D,F,color ="red",title = name,xlabel = "disp",ylabel = "force"))
    sleep(1)

    using Dierckx
    sp1 = Spline1D(D,F)

end
