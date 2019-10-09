// Iscander.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include "Runge.h"
#include "status.h"
//TODO: Эйлер
//TODO: вывод в файл
int main()
{

	status Rocket;

	//status* rocket = new status();
	Rocket.nonIntegr();

	std::vector <status> traj;
	
		
	
	
	status buf;
		while (Rocket.getParam()[4] > 0) {
			traj.push_back(Rocket);
			buf = Rocket;
			runge(Rocket, H);		
		}
		Rocket = buf;

		double h1 = H * 0.5;

		while (abs(Rocket.getParam()[4]) > 1E-3) {
			buf = Rocket;		
			runge(Rocket, h1);
			if (Rocket.getParam()[4] < 0) {
				Rocket = buf;
				h1 /= 2;
			}
		}
		traj.push_back(Rocket);
	std::string afilename = "res.txt";
	std::ofstream fout1(afilename);
	for (int i = 0; i < traj.size(); i++) {
		if ((fabs(traj[i].getParam()[13] *100 - round(traj[i].getParam()[13]*100)) < EPS1) 
			&& ((int)round(traj[i].getParam()[13] * 100) % 2 ==0) 
			|| (i == traj.size()-1)) {
		traj[i].printParam(fout1);
		}
	}
	//double t = H;
	//while ((t < T_FIN) || (abs(t - T_FIN) < EPS1)) {
	//	euler(Rocket, H);
	//	//std::cout << Rocket.getParam()[4] << '\n';
	//	std::cout << Rocket.getParam()[4] << '\t' << Rocket.getParam()[3] << '\n';
	//	Rocket.printParam(fout1);
	//	traj.push_back(Rocket);
	//	t += H;
	//}



		

	
    std::cout << "Hello World!\n"; 
	std::cout << traj[traj.size() - 1].getParam()[5];
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
