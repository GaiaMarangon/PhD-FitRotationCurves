function [rhoNames, factor4pi, nameFactor4pi] = get_modelParams(nameSource,normFactor)

    % select sources
    switch nameSource
        case "Exp"
            rhoNames = ["Exp";"Exp";"Exp"];
        case "TruncatedPlummer"
            rhoNames = ["TruncatedPlummer";"TruncatedPlummer";"TruncatedPlummer"];
        case "HardBall"
            rhoNames = ["HardBall";"HardBall";"HardBall"];
        otherwise
            fprintf("Error: source name not defined\n")
            rhoNames="";
    end

    % select model
    switch normFactor
        case "unitary"
            factor4pi = 1;
            nameFactor4pi = "unitary";
        case "4pi"
            factor4pi = 4*pi;
            nameFactor4pi = "4pi";
        otherwise
            fprintf("Error: model name not defined\n")
            factor4pi = 0;
            nameFactor4pi = "";
    end

end