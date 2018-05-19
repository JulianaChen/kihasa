for x = 1:1:20

    for row_x1 = 1:1:90

        A = S.SS_A(row_x1)
        H = S.SS_H(row_x1)
        K = S.SS_K(row_x1)

        tic
        Base=kron(chebpoly_base(S.nA+1, S.d_A*(A - S.extmin_A) - 1),kron(chebpoly_base(S.nH+1, S.d_H*(H - S.extmin_H) - 1),chebpoly_base(S.nK+1, S.d_K*(K - S.extmin_K) - 1)));
        TVF_poly(row_x1) = sum(coeff(x,:).*Base,2); %cheby_approx
        toc

        tic
        TVF_linear(row_x1) = interpn(unique(S.SS_K),A_wide,unique(S.SS_H),Emax_rsp(:,:,:,x),K,A,H);
        toc
    end

compare(:,:,x) = [TVF(:,x),TVF_poly',TVF_linear'];

end
x = 1:M^3;
plot(x,ABC_func)
plot(x,ABC_func,x,linear)
plot(x,ABC_func,x,linear,x,EV)
plot(x,ABC_func,x,linear,x,EV,x,EV2)