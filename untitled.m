for shuffle=0:2
    for MonkeyNum=1:2
        pSignificant_1_Comb{shuffle+1,MonkeyNum}=pValueSepToComb(pSignificant_Theo1_F{shuffle+1,MonkeyNum},pSignificant_Sim1_F{shuffle+1,MonkeyNum});
        pSignificant_sq_Comb{shuffle+1,MonkeyNum}=pValueSepToComb(pSignificant_Theo2sq_F{shuffle+1,MonkeyNum},pSignificant_Sim2sq_F{shuffle+1,MonkeyNum});
        pSignificant_cross_Comb{shuffle+1,MonkeyNum}=pValueSepToComb(pSignificant_Theo2cr_F{shuffle+1,MonkeyNum},pSignificant_Sim2cr_F{shuffle+1,MonkeyNum});
    end
end

function pComb=pValueSepToComb(px,py)
    pComb=0.5.*erfc(-sqrt((erfcinv(2.*px)).^2+(erfcinv(2.*py)).^2));
end
