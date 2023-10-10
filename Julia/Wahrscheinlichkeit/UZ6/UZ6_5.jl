using Random, Distributions, CairoMakie, StatsBase
#Intialisieren der Variablen und insbesonderen der Speichertypen
#N - Anzahl der Walker
#x -  Vektor der N Mal die 0 enthält
#i - globale Laufvariable
#y - Vektor in dem Anzahl der Versuche bis 0 Punkt erreicht wird
#abgespeichert wird
N = Int64(10000);
global x = zeros(N);
global i = Int64(0);
global y = Float16[];
imax = 100;
#Random walk - gleiche Wahrsch. um in beide Richtungen zu gehen
while i < imax && size(x) != 0
    global i += 1
    global x = x + rand((-1, 1), length(x))
    # Falls die 0 im Vektor x enthalten Intialisieren
    if 0 in x
        #Index der Nullen
        local indic = findall(iszero, x)
        #Initialisierung eine Vektors time der mit Hilfe 
        # der globalen Variable i die Anzahl der Schritte für
        # x Walker speichert
        local time = zeros(length(x))
        local time[indic] .= i
        #Vektor y an dem nach jeder Iteration der Schleife
        # alle Einträge von y die nicht 0 sind gehängt werden
        global y = vcat(y, filter(!iszero, time))
        # Entfernen der Walker die schon einmal 0 erreicht haben
        x = filter(!iszero, x)
    end

end
#Laut Aufgabenstellung wird Mittelwert gebraucht?
mean_y = mean(y);

#Berechnung der Charakteristischen Funktion
#Leider irgendwie echter pain in Julia
myBinomialPDF = N * (binomial.(2 * BigInt.(1:1:imax),BigInt.(1:1:imax)).*(1/4).^BigInt.(1:1:imax) ./ (2 .* BigInt.(1:1:imax).-1))
#Unique Values(Keys) of y mit deren Vielfachheit(Values)
ValuesOfReturn = values(sort(countmap(y)));
KeysOfReturn = keys(sort(countmap(y)));



#PLOT
f = Figure()
Axis(f[1,1],xlabel="x - Steps", ylabel="n zurückgekommen bei x - Steps")
Axis(f[1,2],xlabel="x - Steps", ylabel="Erstes Mal zurückgekommen nach x - Steps")
#Plot CDF
barplot!(f[1,1],Int.(KeysOfReturn),cumsum(ValuesOfReturn);)
lines!(f[1,1],cumsum(myBinomialPDF);color="red")

#Plot Hist
hist!(f[1,2],y;bins = imax, normalization=:none)
lines!(f[1,2],myBinomialPDF;color="red")
f