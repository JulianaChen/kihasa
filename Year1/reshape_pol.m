function [c_func_rsp,m_func_rsp,lr_func_rsp,ln_func_rsp,lu_func_rsp]= reshape_pol(c_func,m_func,lr_func,ln_func,lu_func]

for t=1:1:G.n_period-1
    for x=1:1:G.n_matstat*G.n_wrkexp
        c_func_rsp(:,:,:,:,:,x,t) = reshape(c_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        c_func_rsp9(:,:,:,:,x,t) = reshape(c_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        m_func_rsp(:,:,:,:,:,x,t) = reshape(m_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        m_func_rsp9(:,:,:,:,x,t) = reshape(m_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        lr_func_rsp(:,:,:,:,:,x,t) = reshape(lr_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        lr_func_rsp9(:,:,:,:,x,t) = reshape(lr_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        ln_func_rsp(:,:,:,:,:,x,t) = reshape(ln_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        ln_func_rsp9(:,:,:,:,x,t) = reshape(ln_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
        lu_func_rsp(:,:,:,:,:,x,t) = reshape(lu_func(:,x,t), [G.n_childK,G.n_assets,G.n_hwages,3,3]);
        lu_func_rsp9(:,:,:,:,x,t) = reshape(lu_func(:,x,t),[G.n_childK,G.n_assets,G.n_hwages,G.n_shocks]);
    end
end

