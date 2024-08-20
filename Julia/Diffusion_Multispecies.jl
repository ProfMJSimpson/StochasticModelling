using Plots, DifferentialEquations, LaTeXStrings
gr()


function diff!(du,u,p,t)
dx,N,D1,D2=p
s = zeros(N)

for i in 1:N
s[i]= u[1,i]+u[2,i]
end

for i in 2:N-1
du[1,i]=D1*((2-s[i]-s[i+1])*(u[1,i+1]-u[1,i])-(2-s[i]-s[i-1])*(u[1,i]-u[1,i-1]))/(2*dx^2)+D1*((u[1,i]+u[1,i+1])*(s[i+1]-s[i])-(u[1,i]+u[1,i-1])*(s[i]-s[i-1]))/(2*dx^2)
du[2,i]=D2*((2-s[i]-s[i+1])*(u[2,i+1]-u[2,i])-(2-s[i]-s[i-1])*(u[2,i]-u[2,i-1]))/(2*dx^2)+D2*((u[2,i]+u[2,i+1])*(s[i+1]-s[i])-(u[2,i]+u[2,i-1])*(s[i]-s[i-1]))/(2*dx^2)
end #Exclusion

#for i in 2:N-1
#du[1,i]=D1*(u[1,i+1]-2*u[1,i]+u[1,i-1])/dx^2
#du[2,i]=D2*(u[2,i+1]-2*u[2,i]+u[2,i-1])/dx^2
#end #Nonexclusion



du[1,1]=0.0
du[1,N]=0.0
du[2,1]=0.0
du[2,N]=0.0
end

    
  
    


function pdesolver(L,dx,N,T,C0,D1,D2)
p=(dx,N,D1,D2)
tspan=(0.0,maximum(T))
prob=ODEProblem(diff!,C0,tspan,p)
sol=solve(prob,saveat=T);
return sol
end





LX=400
T=[0,100,200]
dx=0.50
N=Int(LX/dx)+1
x=LinRange(-LX/2,LX/2,N)
D1=1.0;
D2=1.0;

A0=zeros(N)
B0=zeros(N)
AN=zeros(length(T),N)
AB=zeros(length(T),N)
sol=zeros(2,N,length(T))

for i in 1:N
    if x[i] >= -50 && x[i] <=  0
    A0[i]=1
    end
    if x[i] > 0 && x[i] <= 50
        B0[i]=1
    end
end

C0=zeros(2,N)
C0[1,:]=A0
C0[2,:]=B0
sol=pdesolver(LX,dx,N,T,C0,D1,D2); 

p1=plot(x,sol[1,:,1],lw=2,xlims=(-LX/2,LX/2),ylims=(0,1.1), ylabel=L"A(x,t), \, B(x,t)",legend=false)
p1=plot!(x,sol[2,:,1],lw=2,xlims=(-LX/2,LX/2),ylims=(0,1.1), color=:red,legend=false)
p1=plot!(xlims=(-LX/2,LX/2),xticks=([-200,-100,0,100,200],[L"-200",L"-100",L"0",L"100",L"200"]))
p1=plot!(ylims=(0,1.1),yticks=([0,0.25,0.5,0.75,1.0],[L"0.00", L"0.25", L"0.50", L"0.75", L"1.00"]))
p1=plot!(xguidefontsize=12, yguidefontsize=12,xtickfontsize=12, ytickfontsize=12)

p2=plot(x,sol[1,:,2],lw=2,xlims=(-LX/2,LX/2),ylims=(0,1.1), ylabel=L"A(x,t), \, B(x,t)",legend=false)
p2=plot!(x,sol[2,:,2],lw=2,xlims=(-LX/2,LX/2),ylims=(0,1.1), color=:red,legend=false)
p2=plot!(xlims=(-LX/2,LX/2),xticks=([-200,-100,0,100,200],[L"-200",L"-100",L"0",L"100",L"200"]))
p2=plot!(ylims=(0,1.1),yticks=([0,0.25,0.5,0.75,1.0],[L"0.00", L"0.25", L"0.50", L"0.75", L"1.00"]))
p2=plot!(xguidefontsize=12, yguidefontsize=12,xtickfontsize=12, ytickfontsize=12)


p3=plot(x,sol[1,:,3],lw=2,xlims=(-LX/2,LX/2),ylims=(0,1.1), ylabel=L"A(x,t), \, B(x,t)",legend=false)
p3=plot!(x,sol[2,:,3],lw=2,xlims=(-LX/2,LX/2),ylims=(0,1.1), xlabel=L"x" ,color=:red,legend=false)
p3=plot!(xlims=(-LX/2,LX/2),xticks=([-200,-100,0,100,200],[L"-200",L"-100",L"0",L"100",L"200"]))
p3=plot!(ylims=(0,1.1),yticks=([0,0.25,0.5,0.75,1.0],[L"0.00", L"0.25", L"0.50", L"0.75", L"1.00"]))
p3=plot!(xguidefontsize=12, yguidefontsize=12,xtickfontsize=12, ytickfontsize=12)
p4=plot(p1,p2,p3,layout=(3,1))
display(p4)


