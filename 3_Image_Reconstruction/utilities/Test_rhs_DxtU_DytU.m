function RHS = Test_rhs_DxtU_DytU(Wx,Wy,bx,by,tau)

    RHS = tau*(DxtU(Wx-bx)+DytU(Wy-by));

    % compute D'_x(U)
    function dxtu = DxtU(U)
        dxtu = [U(:,end)-U(:, 1) U(:,1:end-1)-U(:,2:end)];
    end

    % compute D'_y(U)
    function dytu = DytU(U)
        dytu = [U(end,:)-U(1, :); U(1:end-1,:)-U(2:end,:)];
    end

end