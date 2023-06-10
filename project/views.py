from django.shortcuts import redirect, render
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO
import base64
from scipy.special import gamma, factorial

# Create your views here.

plot_buffer_list = []
plot_buffer_dict = {}

def index(request):

    global plot_buffer_list, plot_buffer_dict

    if request.method == 'POST':

        C0 = float(request.POST.get('C0', 0.1))
        vm = float(request.POST.get('vm', 1e-4))
        Dm = float(request.POST.get('Dm', 1e-5))
        n = int(request.POST.get('n', 20))
        tau = float(request.POST.get('tau', 1))
        h = float(request.POST.get('h', 0.1))
        tmax = int(request.POST.get('tmax', 1600))
        w = float(request.POST.get('w', 1e-5))  # qarash kerak

        tetim = float(request.POST.get('tetim', 0.1))
        tetm = float(request.POST.get('tetm', 0.1))
        gam = float(request.POST.get('gam', 0.1))

        alpha = float(request.POST.get('alpha', 1))
        bet = float(request.POST.get('bet', 0.1))

        A = 1 + (w*gamma(2-alpha)*tau**alpha)/(gam*tetim)

        Cm = np.zeros((int(tmax) + 1, n+1))
        Cim = np.zeros((int(tmax) + 1, n+1))

        #Chegaraviy shartlar
        for k in range(tmax + 1):
            Cm[k, 0] = C0

        # asosiy qism
        for k in range(tmax):
            for i in range(n+1):
                s1=0
                for l1 in range (k):
                    if k>0:
                        s1+=(Cim[l1+1,i]-Cim[l1,i])*((k-l1+1)**(1-alpha)-(k-l1)**(1-alpha))                       
                Cim[k+1,i]=(Cim[k,i]-s1+(gamma(2-alpha)*tau**alpha*w)*Cm[k,i]/(gam*tetim))/A

            #for i in range(n):
            # Cim[k + 1, i] = A1*(w*(Cm[k,i]-Cim[k,i]))+alpha*Cim[k,i]
                
            for i in range(1, n):
                if bet == 2:
                    s01 = Cm[k, i + 1] - 2 * Cm[k, i] + Cm[k, i - 1]
                else:
                    s01 = 0
                    for l in range(i):
                        s01  = s01 + ((l+1) ** (2 - bet) - (l) ** (2 - bet)) * (Cm[k, i + 1 - l] - 2 * Cm[k, i - l] + Cm[k, i - 1 - l])
                Cm[k + 1, i] = (tau*Dm*s01)/(gamma(3-bet)*(h**bet))-tau*vm*(Cm[k,i]-Cm[k,i-1])/(h)-gam*tau*tetim*s1/(tetm*gamma(2-alpha)*tau**alpha)*(Cim[k+1,i]-Cim[k,i])+Cm[k,i]
            Cm[k + 1, n] = Cm[k + 1,n-1]

        # natijani chop etish
        x = [i * h for i in range(n + 1)]
        for k in range(tmax + 1):
            Cm1 = np.zeros(n+1)
            Cim1 = np.zeros(n+1)
            for i in range(n + 1):
                # print(Cm[k, i], end='  ')
                # print(Cim[k, i], end='  ')
                Cm1[i] = Cm[k, i]
                Cim1[i] = Cim[k, i]
            # print()
            if k % tmax == 0 and k!=0:
                plt.plot(x, Cm1)
                plt.plot(x, Cim1)

                # Templatega chiqarish uchun bufferga saqlanayabdi
                buffer = BytesIO()
                plt.savefig(buffer, format='png')
                buffer.seek(0)
            
        plot_data = base64.b64encode(buffer.read()).decode('utf-8')

        plot_dict = {
            'C0': C0,
            'vm': vm,
            'Dm': Dm,
            'n': n,
            'tau': tau,
            'h': h,
            'tmax': tmax,
            'w': w,
            'tetim': tetim,
            'tetm': tetm,
            'gam': gam,
            'alpha': alpha,
            'bet': bet,
            'plot_data': plot_data
        }

        plot_buffer_list.append(plot_dict)

        context = {
            "C0": C0,
            "vm": vm,
            "Dm": Dm,
            "n": n,
            "tau": tau,
            "h": h,
            "tmax": tmax,
            "w": w,
            "tetim": tetim,
            "tetm": tetm,
            "gam": gam,
            "alpha": alpha,
            "bet": bet,
            'plot_data': plot_data,
            'plot_buffer_list': plot_buffer_list
        }

        return render(request, 'index.html', context)
    return render(request, 'index.html')


def clear(request):
    global plot_buffer_list, plot_buffer_dict
    plot_buffer_list.clear()
    plot_buffer_dict.clear()
    return redirect('index')
