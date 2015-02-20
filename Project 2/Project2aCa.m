% Joel Lubinitsky
% AEE 342 - Project 2a: Analysis of Non-symmetric Airfoil Flows
% Integrate Cp for Ca
% 02/20/15

function area = Project2aCa(n_pan, xCpUpper, xCpLower, CpUpper, CpLower, yUpper, yLower)

    areaUpper = 0;
    areaLower = 0;
    
    for n = [1 : (n_pan / 2) - 2]
        
        dxUpper   = xCpUpper(n + 1) - xCpUpper(n);
        dydxUpper = (yUpper(n + 1) - yUpper(n)) / dxUpper;
        rectUpper = dydxUpper * CpUpper(n) * dxUpper;
        triUpper  = 0.5 * dydxUpper * (CpUpper(n + 1) - CpUpper(n)) * dxUpper;
        areaUpper = areaUpper + rectUpper + triUpper;
        
        dxLower   = xCpLower(n + 1) - xCpLower(n);
        dydxLower = (yLower(n + 1) - yLower(n)) / dxLower;
        rectLower = dydxLower * CpLower(n) * dxLower;
        triLower  = 0.5 * dydxLower * (CpLower(n + 1) - CpLower(n)) * dxLower;
        areaLower = areaLower + rectLower + triLower;
        
    end
    
    area = areaUpper - areaLower;