function [nameGal,indexGal] = get_galaxyParams(nameGal)

    % select galaxy index for experimental data loading
    switch nameGal
        case "UGC02953"
            indexGal = 110;
        case "NGC5055"
            indexGal = 82;
        case "UGC09037"
            indexGal = 162;
        otherwise
            fprintf("Error: galaxy name not defined\n")
            indexGal = 0;
    end
end