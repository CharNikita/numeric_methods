#include "MatrixSystem.cpp"

#include <fstream>
#include <iomanip>
#include <ctime>

using namespace std;

int main()
{
   //MatrixSystem<double> slv = MatrixSystem<double>("params", "ia_gil"); //Гильберт
   //MatrixSystem<double> slv = MatrixSystem<double>("param_b2", "ia_b2", "ja_b2", "di_b2", "au_b2", "al_b2", "f_b2");
   //MatrixSystem<double> slv = MatrixSystem<double>("params1", "ia1", "ja1", "di1", "au1", "al1","f1");
   //MatrixSystem<double> slv = MatrixSystem<double>("param_b1", "ia_b1", "ja_b1", "di_b1", "au_b1", "al_b1", "f_b1");

   //MatrixSystem<double> gil = MatrixSystem<double>("kuslau", "gi", 1);
   //return 0;

   bool is_done = true;
   int flag;

   while (is_done)
   {
      cout << "Methods" << endl;
      cout << "1. MSG" << endl;
      cout << "2. LOS" << endl;
      cout << "3. MSG_sq" << endl;
      cout << "4. LOS_sq" << endl;
      cout << "5. MSG_di" << endl;
      cout << "6. LOS_di" << endl;
      cout << "Enter the number of method:" << endl;
      cin >> flag;

      if (flag >= 1 && flag <= 6) is_done = false;
      else cout << "Enter the correct number." << endl;
   }

   MatrixSystem<double> ms = MatrixSystem<double>("kuslau",
      "gi", "gj", "di", "ggu", "ggl", "pr", flag);
   cout << setprecision(15) << scientific;

   clock_t start = clock();

   switch (flag)
   {
      case(1):
      {
         cout << "MSG" << endl;
         ms.msg_basic();
         break;
      }
      case(2):
      {
         cout << "LOS" << endl;
         ms.los_basic();
         break;
      }
      case(3):
      {
         cout << "MSG_sq" << endl;
         ms.msg_sq();
         break;
      }
      case(4):
      {
         cout << "LOS_sq" << endl;
         ms.los_sq();
         break;
      }
      case(5):
      {
         cout << "MSG_di" << endl;
         ms.msg_di();
         break;
      }
      case(6):
      {
         cout << "LOS_di" << endl;
         ms.los_di();
         break;
      }
      default: break;
   }

   clock_t finish = clock();

   ofstream fout;
   fout.open("output.txt");
   fout << scientific << setprecision(15);
   for (int i = 0; i < ms.getSize(); i++)
      fout << ms.x[i] << endl;
   fout.close();

   clock_t time = finish - start;
   cout << "time=" << time << endl;
   //cout << "nev=" << << endl;

   return 0;
}
