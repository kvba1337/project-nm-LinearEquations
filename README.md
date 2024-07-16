Metody Numeryczne

Projekt 2 – Układy równań liniowych

dr hab. inż. Grzegorz Fotyga, prof. PG

24 marca 2024

1. **Wstęp**

Celem projektu jest implementacja i analiza dwóch metod iteracyjnych (Jacobiego i Gaussa-Seidla) oraz jednej metody bezpośredniej (faktoryzacja LU) rozwiązywania układów równań liniowych. Testy poszczególnych metod będą przeprowadzane na układach równań, które mogą powstać w wyniku dyskretyzacji równań różniczkowych i są powszechnie stosowane w takich zagadnieniach jak: elektronika, elektrodynamika, mechanika (zastosowania lotnicze, biomechanika, motoryzacja), badanie wytrzymałości materiałów i konstrukcji, symulacje odkształceń, naprężeń, przemieszczeń

i drgań, akustyka, fotonika, termodynamika, dynamika płynów i wiele innych.

W rzeczywistych problemach rozwiązywane są układy równań zawierające setki milionów nie- wiadomych, dla których obliczenia trwają często wiele godzin, a nawet dni, mimo wykorzystywania najnowszych superkomputerów. Opracowanie nowych efektywnych metod rozwiązań (dostosowanych do współczesnych architektur komputerowych) jest dużym wyzwaniem zarówno z punktu widzenia matematyki, jak i informatyki. Jest ono przedmiotem badań wielu ośrodków naukowych, ponieważ bez niego rozwój wymienionych wyżej dziedzin wiedzy byłby **niemożliwy** .

W praktyce najczęściej stosuje się tak zwany rzadki format przechowywania macierzy (który przechowuje tylko wartości niezerowe i ich położenie w macierzy), ponieważ zdecydowana większość elementów ma wartość 0. Jednak ze względu na prostotę, w ramach tego projektu domyślnie będzie wykorzystywany tzw. format pełny (przechowujący wszystkie wartości, również 0), który może być stosowany w problemach zawierających zazwyczaj nie więcej niż kilka tysięcy niewiadomych. Mimo, że testowane będą jedynie podstawowe metody rozwiązań, wykonanie projektu będzie dobrym fun- damentem do poznania bardziej zaawansowanych metod iteracyjnych (np. metody gradientów sprzę- żonych, GMRES, QMR itp.).

2. **Konstrukcja układu równań**

Układ równań liniowych ma następującą postać:

<a name="_page0_x278.48_y597.44"></a>**Ax** = **b** (1) gdzie **A** jest macierzą systemową [^1][,](#_page0_x61.59_y741.36) **b** jest wektorem pobudzenia [^2][,](#_page0_x61.59_y765.27) natomiast **x** jest wektorem rozwią- zań reprezentującym szukaną wielkość fizyczną [^3][.](#_page0_x61.59_y777.22)

- Na potrzeby projektu należy przyjąć, że **A** jest tzw. macierzą pasmową o rozmiarze *N* × *N* i zdefiniowaną w ([2),](#_page1_x166.60_y71.15) gdzie *N* ma wartość 9*cd*, *c* jest przedostatnią cyfrą numeru Twojego indeksu, natomiast *d* ostatnią (np. dla indeksu 102263 *N* = 963). Macierz **A** zawiera więc pięć diagonali - główna z elementami *a*1, dwie sąsiednie z elementami *a*2 i dwie skrajne diagonale z elementami *a*3.

2



*a*1 

 *a*2 

<a name="_page1_x166.60_y71.15"></a>**A** =  *a*03

 ...

*a*2 *a*3 0 *a*1 *a*2 *a*3 *a*2 *a*1 *a*2

*a*3 *a*2 *a*1 ... ... ...



0 0 0 *...* 0

0 0 0 *...* 0 

*a*3 0 0 *...* 0 

*a*2 *a*3 0 *...* 0 *,* (2)



.. .. .. .. .. 

. . . . . 





0 0 *...* 0 0 0 *a*3 *a*2 *a*1

- Prawa strona równania to wektor **b** o długości *N* .
- W wyniku rozwiązania układu równań (1[) otrz](#_page0_x278.48_y597.44)ymujemy wektor **x**.
3. **Wektor residuum**

Ważnym elementem algorytmów iteracyjnych (np. Jacobiego i Gaussa-Seidla) jest określenie w której iteracji algorytm powinien się zatrzymać. W tym celu najczęściej korzysta się z residuum [1], czyli wektora który dla *k* – tej iteracji przyjmuje postać:

**res**(*k*) = **Ax** (*k*) − **b***.* (3) Badając normę euklidesową residuum ( *norm* (**res**(*k*))), możemy w każdej iteracji algorytmu obli-

czyć jaki błąd wnosi wektor **x**(*k*). Jeżeli algorytm zbiegnie się do dokładnego rozwiązania, to residuum stanowić będzie wektor zerowy. Ponieważ w praktyce osiągnięcie dokładnego rozwiązania metodami iteracyjnymi jest niespotykane, to jako kryterium stopu przyjmuje się osiągnięcie normy residuum

mniejszej niż np. 10− 6.

Residuum nazywane jest również wektorem reszt [2[\] lub](#_page2_x43.65_y73.64) wektorem residualnym [3]. W literaturze anglojęzycznej residuum określane jest jako *residual vector* lub *residual* [4].

4. **Zadania**

Sprawozdanie powinno zawierać m.in. analizę rezultatów osiągniętych w zadaniach **B**, **C**, **D**, **E**.

- **Zadanie A** – Stwórz układ równań dla *a*1 = 5+ *e*, gdzie *e* jest czwartą cyfrą Twojego indeksu, *a*2 = *a*3 = −1. Rozmiar macierzy *N* zdefiniowano w punkcie 2 tej instrukcji. **b** jest wektorem
  - długości *N* , którego *n*−ty element ma wartość *sin* (*n* ·(*f* + 1)) , gdzie *f* jest trzecią cyfrą Twojego indeksu. We wstępie sprawozdania opisz rozwiązywane równanie macierzowe. (5%)
- **Zadanie B** – Zaimplementuj metody iteracyjne rozwiązywania układów równań liniowych: Jacobiego i Gaussa–Seidla. Opisz ile iteracji potrzebuje każda z nich w celu wyznaczenia roz- wiązania układu równań z zadania **A**, którego norma residuum jest mniejsza niż 10− 9. Dla obu metod przedstaw na **wykresie** jak zmienia się norma residuum w kolejnych iteracjach wykony- wanych w celu wyznaczenia rozwiązania (oś *y* w skali logarytmicznej). Porównaj czas trwania algorytmów. (30%)
- **Zadanie C** – Stwórz układ równań dla *a*1 = 3, *a*2 = *a*3 = −1, natomiast *N* i wektor **b** określ zgodnie z treścią zadania **A**. Czy metody iteracyjne dla takich wartości elementów macierzy

  **A** zbiegają się? Dla obu metod przedstaw na **wykresie** jak zmienia się norma residuum w kolejnych iteracjach (oś *y* w skali logarytmicznej). (10%)

- **Zadanie D** – Zaimplementuj metodę bezpośredniego rozwiązania układów równań liniowych: metodę faktoryzacji LU i zastosuj do równania badanego w p. **C**. Ile wynosi norma residuum w tym przypadku? (30%)
- **Zadanie E** – Stwórz **wykres** zależności czasu wyznaczenia rozwiązania dla trzech badanych metod w zależności od liczby niewiadomych *N* = {100*,*500*,*1000*,*2000*,*3000*...* } dla macierzy opisanej w zadaniu **A**. (10%)
- **Zadanie F** – Zwięźle opisz swoje obserwacje po wykonaniu zadań **A**–**E**. (15%)

**Literatura**

1. Bjorck<a name="_page2_x43.65_y55.70"></a> A., Dahlquist G., *Metody numeryczne* , PWN, 1987
1. Fortuna<a name="_page2_x43.65_y73.64"></a> Z., Macukow B., Wąsowski J., *Metody numeryczne* , PWN, 2017
1. Kincaid,<a name="_page2_x43.65_y93.06"></a> Cheney, *Analiza numeryczna* , WNT, 2006
1. Saad,<a name="_page2_x43.65_y112.49"></a> Yousef, *Iterative Methods for Sparse Linear Systems* , Society for Industrial and Applied Mathematics, 2003.
4

[^1]: <a name="_page0_x61.59_y741.36"></a>W zależności od problemu może ona reprezentować np. obwód elektroniczny, geometrię sali koncertowej, turbinę, karoserię samochodu itp.
[^2]: <a name="_page0_x61.59_y765.27"></a>np. impuls elektroniczny, wektor siły, fala dźwiękowa itp.
[^3]: <a name="_page0_x61.59_y777.22"></a>np. rozkład pola elektromagnetycznego, natężenie dźwięku itp.