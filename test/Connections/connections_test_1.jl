#Put your code here Sangchu.  
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
el = DataFrame(CSV.File("D:/JHU Spring 2021/Fastener Testing Data_Zhidong/Trail1/Mon.csv"))

# Convert the already read CSV file into array
df_matrix = Matrix(el)

# Visualize the chosen node
for i in 1:size(df_matrix,1)

    # Load the 4 node's displacement and force
    F = df_matrix[i,6:9]
    D = df_matrix[i,2:5]
    #
    name = df_matrix[i,1]

    df = DataFrame(force = F, displacement = D)
    plot(D,F)
    Plots.display(plot(D,F,color ="red",title = name,xlabel = "disp",ylabel = "force"))
    sleep(1)
    println(i)
end
