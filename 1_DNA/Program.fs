﻿// Learn more about F# at http://fsharp.net
// See the 'F# Tutorial' project for more help.
open System
open System.IO

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)

type FASTA =
    {Label:string;String:string}
    
let parseFasta (text:string) =
    text.Split([|'>'|],StringSplitOptions.RemoveEmptyEntries)
    |> Seq.map 
        (fun lines-> 
            lines.Split([|'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries)
            |> fun toks -> 
                {Label=toks.[0]; String=toks.[1..]|>String.concat ""})

let guard m c =
    if Map.containsKey c m then
        m.[c]
    else
        0


// 1. DNA

getData "dna"
|> Seq.countBy id
|> Map.ofSeq
|> fun m -> printfn "%d %d %d %d" (guard m 'A') (guard m 'C') (guard m 'G') (guard m 'T')

// 2. RNA: Transcribing DNA into RNA
let dnaToRna (x:string) : string = x.Replace('T','U')

getData "rna"
|> dnaToRna
|> printfn "%s"


// 3. REVC: reverse complement
let complement c =
    match c with
    | 'A' -> 'T'
    | 'T' -> 'A'
    | 'C' -> 'G'
    | 'G' -> 'C'
    | _ -> raise (new Exception("Invalid nucleotide"))

let reverseComplement x =
    x
    |> Seq.map complement 
    |> List.ofSeq
    |> List.rev
    |> String.Concat

getData "revc"
|> fun x->x.Trim()
|> reverseComplement
|> printfn "%s"

// 4. FIB: fibonnaci rabbits
let rec pairsAfterN n k =
    if n=1 || n=2 then 1L
    else 
        let current=pairsAfterN (n-1) k
        let offspring=k*(pairsAfterN (n-2) k)
        current+offspring

getData "fib"
|> fun l -> l.Split(' ')
|> fun toks -> pairsAfterN (Int32.Parse(toks.[0])) (Int64.Parse(toks.[1]))
|> printfn "%d"


// 5. FIBD: mortal fib rabbits
// todo: handle M
let rec pairsAfterNwithM n m =
    if n=1 then [1L]
    else if n=2 then [0L;1L]
    else
        let lastMonth=
            pairsAfterNwithM (n-1) m            
        (lastMonth|>Seq.skip 1|>Seq.sum) :: (lastMonth|>Seq.truncate (m-1)|>List.ofSeq)

getData "fibd"
|> fun l -> l.Split(' ')
|> fun toks -> pairsAfterNwithM (Int32.Parse(toks.[0])) (Int32.Parse(toks.[1]))
|> Seq.sum
|> printfn "%d"


// 6. GC: computing GC content

let gcContent s =
    let m = Seq.countBy (fun x -> x) s |> Map.ofSeq
    100.*(((guard m 'C')+(guard m 'G'))|>float)/(Seq.length s|>float)


getData "gc"
|> fun x -> x.Trim()
|> parseFasta
|> Seq.maxBy (fun f -> gcContent f.String)
|> fun f -> printfn "%s\n%f" f.Label (gcContent f.String)


// 7. HAMM: counting point mutations (hamming distance)
let dh s1 s2 =
    Seq.zip s1 s2
    |> Seq.map (fun (c1,c2) -> if c1=c2 then 0 else 1)
    |> Seq.sum


getData "hamm"
|> fun x -> x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)
|> fun lines -> printfn "%d" (dh lines.[0] lines.[1])

// 8. IPRB: Mendel's first law.
let rec partFact (n:int64) (stop:int64) =
    if (n<=stop) then 1L
    else n*(partFact (n-1L) stop)

let rec fact (n:int64) = partFact n 1L 
//let C n k = (partFact n (n-k))/(fact k) // blows up for k>=67
// https://en.wikipedia.org/wiki/Binomial_coefficient#Recursive_formula
// too expensive.
(*let rec C n k = 
    if k=0L || k=n then 1L
    else (C (n-1L) (k-1L)) + (C (n-1L) k)*)
// https://en.wikipedia.org/wiki/Binomial_coefficient#Multiplicative_formula
// This was actually incorrect on the wikipedia page and I had to make an edit!
let rec C n k =
    if k=0L || k=n then 1L
    else if k>=(n/2L+1L) then C n (n-k)
    else 
        [1L..k]
        |> Seq.map (fun i -> (float (n-k+i))/(float i))
        |> Seq.reduce (*)
        |> int64

let pDom k m n =
    let kk = float (C k 2L)
    let km = float (k*m)
    let mm= 0.75*float (C m 2L)
    let kn= float (k*n)
    let mn=0.5*float (m*n)
    (kk+km+mm+kn+mn)/float (C (k+m+n) 2L)
// k: num  homozygous dominant
// m: num heterozygous
// n: num homozygous recessive

let p64 x = Int64.Parse(x)

getData "iprb"
|> fun l -> l.Split(' ')
|> fun t -> pDom (p64 t.[0]) (p64 t.[1]) (p64 t.[2])
|> printfn "%f"

// 9. PROT Translating RNA into Protein
// be sure source is materialized or you'll be re-traversing the sequence over and over...
let chop segmentSize source =
    source    
    |> Seq.unfold(fun chopSource ->
        if Seq.isEmpty chopSource then None else
        let segment = chopSource |> Seq.truncate segmentSize
        let rest = chopSource |> Seq.skip (Seq.length segment)
        Some(segment, rest)
    )

let chopStr segmentsize (source:string) =
    seq {
        for i in [0..segmentsize..source.Length-1] do
            yield source.Substring(i,segmentsize)
    }

let decoder =
    @"UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G".Split([|' ';'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries)
    |> chop 2
    |> Seq.map (fun s->Array.ofSeq s |> fun x->x.[0],x.[1])
    |> Map.ofSeq
           
let decode (s:string) =
    chopStr 3 s
    |> Seq.map (fun x -> x|> Seq.map string|>String.concat "")    
    |> Seq.map (fun x -> decoder.[x])
    |> Seq.takeWhile (fun x -> x<>"Stop")
    |> String.concat ""
            
getData "prot" |> decode|>printfn "%s" 

// 10. SUBS. Finding a Motif in DNA
let allOccurrences s t =
    let rec allOccurrencesHelper (s:String) (t:String) =
        let ix=s.IndexOf(t)
        if (ix>=0) then 
            (ix)::(allOccurrencesHelper (s.Substring(ix+1)) t)
        else
            []
    allOccurrencesHelper s t
    |> Seq.scan (fun s t -> s+t+1) -1
    |> Seq.map ((+) 1)
    |> Seq.skip 1
    |> List.ofSeq
    // test: allOccurrences "GATATATGCATATACTT" "ATAT";;
    // allOccurrences "CTGAGTGGCTGAGTGCTGAGTGCTGAGTGCGAGCTGAGTGCCTGAGTGCTGAGTGCTGAGTGGTAATCCCTGAGTGTCGAACTGAGTGATACTGAGTGACATCTGAGTGCTGAGTGAACCTCTGAGTGCTGAGTGACCTGAGTGTTCTGTCTGAGTGCACGAACAGCTTCCTGAGTGCTTCGCTGAGTGCCTGAGTGATCTGAGTGCTGAGTGCCTGAGTGCTGAGTGTATTGCTGAGTGAAACTGAGTGTCCTGAGTGCTGAGTGCTGAGTGCCTCACTGAGTGTCTCTACTGAGTGCGCTGAGTGGCTGAGTGTGCCTGAGTGTGCTGAGTGTGCTGAGTGTCTGAGTGCTTCCTGAGTGAGGTGACGGCGCTGAGTGCTGAGTGTGCTACTGAGTGCTGAGTGCTGAGTGCATCGCTGAGTGTCCTGAGTGCTGAGTGTCAGCTGAGTGAGAACTGAGTGTCTGAGTGCGCCGCTGAGTGTCTGAGTGAGGAACTGAGTGGGCCTGAGTGAACCTGAGTGCTGAGTGTAAGCTGAGTGACTGAGTGCACTGAGTGACCGCTGAGTGACCTGAGTGCCTGAGTGCTTTAAACTGAGTGTGTGGCGACTGAGTGCTGAGTGGACCTGAGTGTCTGAGTGGCTGAGTGGTGCTATTCTGAGTGCTGAGTGCTGAGTGATCTGAGTGGCCCTGAGTGTCTGAGTGTGTCTAATACTCTGAGTGTTAGCTGAGTGTTTAAGCTGAGTGCTGAGTGCTGAGTGCTGAGTGATGCTGAGTGGCGGGCTGAGTGATCTGAGTGATACTGAGTGCTGAGTGCGCTGAGTGGCTGAGTGCATGCTGAGTGCTGAGTGAGACGACCGAACTGAGTGCTGAGTGACTGAGTGAGTAGACTGAGTG" "CTGAGTGCT" |> Seq.map string |>String.concat " "
getData "subs"
|> fun x->x.Split([|'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries)
|> fun t -> allOccurrences t.[0] t.[1] 
|> Seq.map string 
|> String.concat " "
|> printfn "%s"

// 11. CONS Consensus and profile
#if INTERACTIVE
#r @"c:\GitHub\rosalind\packages\MathNet.Numerics.3.8.0\lib\net40\MathNet.Numerics.dll"
#r @"c:\GitHub\rosalind\packages\MathNet.Numerics.FSharp.3.8.0\lib\net40\MathNet.Numerics.FSharp.dll"
#endif
open MathNet.Numerics.LinearAlgebra

type ProfileMatrix = 
    Map<char,MathNet.Numerics.LinearAlgebra.Vector<double>>
let addProfiles p1 p2 = 
    let keys=
        [Map.toSeq p1 ;Map.toSeq p2]
        |> Seq.concat
        |> Seq.map (fun (k,v) ->k)
        |> Set.ofSeq
    keys
    |> Seq.fold 
        (fun s k -> 
            if Map.containsKey k p1 then
                if Map.containsKey k p2 then
                    s|>Map.add k (p1.[k]+p2.[k])
                else
                    s|>Map.add k p1.[k]
            else
                s|>Map.add k p2.[k])
            Map.empty

let consensus (p:ProfileMatrix) =
    let len = p|>Map.toSeq |> Seq.take 1 |> Seq.exactlyOne |> fun kvp ->p.[fst kvp].Count
    [0..len-1]
    |> Seq.map (fun i -> 
                ['A';'C';'G';'T']
                |> Seq.maxBy (fun c -> p.[c].[i])
                )
    |> Seq.map string
    |> String.concat ""

let printProfVector (v:Vector<double>) =
    v.ToArray()
    |> Seq.map string
    |> String.concat " "

let printProfMtx (p:ProfileMatrix) =
    ['A';'C';'G';'T']
    |> Seq.iter (fun c->printfn "%c: %s" c (printProfVector p.[c]))

let V = Vector<double>.Build
let oneAt len i = V.Dense(len, fun _i -> if _i=i then 1. else 0.)
let makeProfile (f:FASTA) =
    let len = f.String.Length
    let initState=
        ['A';'C';'G';'T']
        |> Seq.map (fun k -> (k,V.Dense(len)))
        |> Map.ofSeq
    [0..len-1]
    |>Seq.fold 
        (fun s i -> 
            let c=f.String.[i]
            Map.empty|>Map.add c (oneAt len i)
            |>addProfiles s
            ) initState


getData "cons"
|> fun x -> x.Trim()
|> parseFasta
|> Seq.map makeProfile
|> Seq.reduce addProfiles
|> fun p -> 
    printfn "%s" (consensus p)
    printProfMtx p

// finds all elements v of V such that v ends with the prefix of f up to length k 
let Ok k V f =
    let str=f.String;
    let prefix=str.Substring(0,k)
    V
    |> Seq.filter (fun v -> v.Label<>f.Label && v.String.EndsWith(prefix))
    |> Seq.map (fun v -> v,f)

//11. GRPH
getData "grph"
|> fun x -> x.Trim()
|> parseFasta
|> fun V -> Seq.map (Ok 3 V) V
|> Seq.concat
|> Seq.iter (fun (v,f)-> printfn "%s %s" v.Label f.Label)

// 12. IEV: expected offspring
// we only need the probs but for troubleshooting include the names
let pTable = 
    [| "AA-AA",1.0; "AA-Aa",1.0; "AA-aa",1.0; "Aa-Aa",0.75; "Aa-aa",0.5; "aa-aa",0.0|]
    |> Array.map (fun (k,v)->v)

let ievdata = 
    getData "iev"
    |> fun x->x.Split([|' '|],StringSplitOptions.RemoveEmptyEntries)
    |> Seq.map (Double.Parse)
    |> Array.ofSeq
ievdata
|> Seq.map2 (*) pTable
|> Seq.sum
|> (*) 2.0
|> printfn "%f"

// 13. LCSM. Finding a shared motif.  
// Wherein impatience gets the better of me and I go for the brute force solution.  Strings are guaranteed to be less than 1kbp.
// My first attempt built sets from every enumerated substring of every string but it was too expensive, so I enumerate just the first string
// and then pare down the set as I check for containment in each of the subsequent strings.
// The 'right' way to do this is with a generalized suffix tree, but the expense of implementing one seemed too high.  Laziness is the programmer's virtue.

let enumerateSubstrings (s:string) =
    seq {
        for i in [0..s.Length-1] do
            for j in [1..s.Length-i] do
                yield s.Substring(i,j)
    }
let substrSet s = enumerateSubstrings s |> Set.ofSeq

let intersectWithString (s:string) (ix:Set<string>) =
    Set.toSeq ix
    |> Seq.filter (fun s' -> s.Contains(s'))
    |> Set.ofSeq

getData "lcsm"
|> fun x -> x.Trim()
|> parseFasta
|> List.ofSeq
|> fun x -> 
    match x with 
    | [] -> []
    | head :: tail ->
        Seq.fold (fun ix s -> intersectWithString s.String ix) (substrSet head.String) tail
        |> Set.toList    
|> List.maxBy (fun x -> x.Length)
|> printfn "%s"

// 14. LIA 
// at least N in k generations.
// This works out to at least N successes in 2^k Bernoulli trials with success probability 0.25.
open System
type Allele = char
type Gene = string
let combineAlleles a1 a2 : string =    
    if Char.ToUpper(a1)=a1 then
        sprintf "%c%c" a1 a2
    else
        sprintf "%c%c" a2 a1

let combineGenes (g1:Gene,p1:float) (g2:Gene,p2:float) : seq<string*float> = // yields the punnet square probabilities for a combination of genes
    seq {
        for a1 in g1 do
            for a2 in g2 do
                yield (combineAlleles a1 a2, p1*p2/4.0)
    }
    |> Seq.groupBy (fun (g,p)->g)
    |> Seq.map (fun (g,c)->g,c|>Seq.sumBy(fun (g',p')->p' ))

let combineGenes2 (gs1:seq<Gene*float>) (gs2:seq<Gene*float>) =
    seq {
        for g1 in gs1 do
            for g2 in gs2 do
                yield combineGenes g1 g2
    }
    |> Seq.concat
    |> Seq.groupBy (fun (g,p)->g)
    |> Seq.map (fun (g,c)->g,c|>Seq.sumBy(fun (g',p')->p' ))

type Genome = Map<string,seq<string*float>> // just a list of genes with probabilities.
let makeGenome (strs:string seq) (ps:float seq) : Genome =
    Seq.zip strs ps
    |> Seq.map (fun (s,p)->Char.ToUpper(s.[0])|>string,[s,p]|>List.toSeq)
    |> Map.ofSeq
    

let mateGenomes (g1:Genome) (g2:Genome) =
    let keys m = Map.toSeq m|> Seq.map (fun (k,v)->k)|> List.ofSeq|>Set.ofList
    let g1keys=keys g1
    let g2keys=keys g2
    let allkeys = Set.union g1keys g2keys
    allkeys
    |> Set.toSeq
    |> Seq.map (fun k -> k,combineGenes2 (g1.[k]) (g2.[k]))
    |> Map.ofSeq

(* tool around with these functions for a bit and see how they interact. 
let g1=makeGenome ["Xx";"Yy"] [1.0;1.0];;

val g1 : Genome = map [("X", [("Xx", 1.0)]); ("Y", [("Yy", 1.0)])]

> mateGenomes g1 g1;;
val it : seq<string * seq<string * float>> =
  seq
    [("X", seq [("XX", 0.25); ("Xx", 0.5); ("xx", 0.25)]);
     ("Y", seq [("YY", 0.25); ("Yy", 0.5); ("yy", 0.25)])]

> let g2=mateGenomes g1 g1;;

val g2 : Map<string,seq<string * float>> = map [("X", <seq>); ("Y", <seq>)]

> mateGenomes g1 g2;;
val it : Map<string,seq<string * float>> =
  map
    [("X", seq [("XX", 0.25); ("Xx", 0.5); ("xx", 0.25)]);
     ("Y", seq [("YY", 0.25); ("Yy", 0.5); ("yy", 0.25)])]
> let g3=mateGenomes g1 g2;;

val g3 : Map<string,seq<string * float>> = map [("X", <seq>); ("Y", <seq>)]

> mateGenomes g3 g1;;
val it : Map<string,seq<string * float>> =
  map
    [("X", seq [("XX", 0.25); ("Xx", 0.5); ("xx", 0.25)]);
     ("Y", seq [("YY", 0.25); ("Yy", 0.5); ("yy", 0.25)])]

The key observation is "Xx" and "Yy" always have 0.5 probability (so, jointly they have 0.25 probability).

Now the problem is just one of the likelihood of producing >=34 heads in 128 bernoulli trials, where p(success)=0.25

*)


// Blows up. Argh      
let binomOld (p:float) n k = 
    let cnk=(C n k)|>float
    printfn "%f %d %d %f" p n k cnk    
    let k' = float k
    let n' = float n
    cnk*(p**k')*((1.-p)**(n'-k'))

let N u s x =
    let coeff=1./(s*sqrt 2.*Math.PI)
    let arg = (x-u)**2./(2.*s**2.)
    coeff * exp (-arg)
//https://en.wikipedia.org/wiki/Binomial_distribution#Normal_approximation
let binom (p:float) (n:int64) (k:int64) =   
    //printfn "%f %f %f" p n k
    if (float n*p>=10. && float n*p*(1.-p)>=10.) then
        N (float n*p) (float n*p*(1.-p)) (float k)
    else
        float (C n k)*(p**(float k))*((1.0-p)**float (n-k))

let probHeteroAtLeastN (p:float) (k:int64) (N:int64) = // p(genome) after k generations. 
    let k' = float k    
    //let p=0.25    
    let n=2.**k'|>int64
    let bnm = binom p n
    [N..(int64 (2.**k'))]
    |> Seq.map (fun i -> i |> bnm)
    |> Seq.sum

getData "lia"
|> fun x -> x.Split(' ')
|> fun t -> probHeteroAtLeastN 0.25 (Int64.Parse(t.[0])) (Int64.Parse(t.[1]))
|> printfn "%f"

// 14. MPRT
// Failed on the first attempt- maybe a copy-paste or download failure?
#if INTERACTIVE
#r @"C:\GitHub\rosalind\packages\Http.fs.1.5.1\lib\net40\HttpClient.dll"
#endif
open HttpClient   
open System.Text.RegularExpressions

let fetchSequence id =  
    createRequest Get (sprintf "http://www.uniprot.org/uniprot/%s.fasta" id) 
    |> getResponseBody
    |> parseFasta
    |> Seq.take 1
    |> (fun s -> id,Seq.exactlyOne s)

// todo: proper regex translation
// Needs help with overlapping regexes. see: http://stackoverflow.com/questions/320448/overlapping-matches-in-regex
let nglycosylation = new Regex("N[^P][S|T][^P]")
let matchLocations (rex:Regex) (str:string) =
    Seq.unfold 
        (fun (ix) -> 
            if ix>=str.Length then // todo: check length.
                None
            else 
                let m = nglycosylation.Match(str,ix) 
                if (m.Success) then
                    Some (m.Index+1,m.Index+1)
                else None) 0
    //let matches=rex.Matches(str)
    //seq { for m in matches do yield (m.Index+1) }


getData "mprt"
|> (fun x-> x.Trim().Split([|'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries))
|> Seq.map fetchSequence
|> Seq.map (fun (id,f)-> id,matchLocations nglycosylation f.String)
|> Seq.iter 
    (fun (f,s)->
        if Seq.isEmpty s then ()
        else
            printfn "%s\n%s" f (s|>Seq.map string|>String.concat " ")
    )

// 15. MRNA
// Explode protein string to count possible RNA strings.  Basically walk along 
let encoder =
    decoder
    |> Map.toSeq
    |> Seq.groupBy (fun (x,y)-> y)
    |> Seq.map (fun (y,s)->(y,Seq.length s))
    |> Map.ofSeq

getData "mrna"
|> fun x->x.Trim()
|> Seq.fold (fun (s:int) (c:char) -> (s*encoder.[string c])%1000000) 1
|> (*) encoder.["Stop"]
|> fun x-> x%1000000
|> printfn "%d"

type ParseState = 
    INIT of string list
    |PARSING of (string list)*(string list)
    
// would unfold work better here?
let parseNext (str:string,i:int,m:Map<int,ParseState>) c =
    //eprintfn "%A" (str,i,m)
    let mapstate=m.[i%3]
    if (i>str.Length-3) then
        (str,i+1,m)
    else
        let workingStr = str.Substring(i,3)
        match mapstate with
        | INIT(l) -> 
            if workingStr="AUG" then
                (str,i+1,m|>Map.add (i%3) (PARSING([decoder.[workingStr]],l)))
            else
                (str,i+1,m)
        | PARSING(s,l) ->
            if decoder.[workingStr]="Stop" then  // found the stop codon
                (str,i+1,m|>Map.add (i%3) (INIT([s;l]|>List.concat)))
            else
                let c=decoder.[workingStr]                    
                let sOut=s|>List.map (fun x->x+c)
                if workingStr="AUG" then // found another start codon. Add another parse string.
                    (str,i+1,m|>Map.add (i%3) (PARSING(c::sOut,l)))
                else
                    (str,i+1,m|>Map.add (i%3) (PARSING(sOut,l)))

let parseProteins (str:string) =
    Seq.fold parseNext (str,0,Map.ofList [0,INIT([]);1,INIT([]);2,INIT([])]) str
    |> fun (_,_,m) -> m
    |> Map.toSeq
    |> Seq.map 
        (fun (ix,p) ->
            match p with
            | INIT(l) -> l
            | PARSING(s,l) -> [l]|>List.concat // was [s;l]
            )
    |> List.concat   
  


getData "orf"
|> fun x->x.Trim()
|> parseFasta
|> Seq.exactlyOne
|> fun x -> x.String
|> fun x->    
    [
        x|>dnaToRna |> parseProteins
        x|> reverseComplement |> dnaToRna |> parseProteins
    ]
|> Seq.concat
|> Seq.groupBy id
|> Seq.map (fun (x,_)->x)
|> Seq.iter (printfn "%s")

// perm
// stolen shamelessly from 
// 1. http://stackoverflow.com/questions/1526046/f-permutations 
// which was copied from p 166-167 of 
// 2. http://www.ffconsultancy.com/products/fsharp_for_technical_computing/?so
// which now I think I should buy...
let rec distribute e = function
  | [] -> [[e]]
  | x::xs' as xs -> (e::xs)::[for xs in distribute e xs' -> x::xs]

let rec permute = function
  | [] -> [[]]
  | e::xs -> List.collect (distribute e) (permute xs)

let perm n =
    printfn "%d" (fact n)
    [1L..n]
    |> permute
    |> Seq.map (fun x-> x|>Seq.map string|>String.concat " ")
    |> Seq.iter (printfn "%s")

// prtm
let monomasstbl = 
    "A   71.03711
C   103.00919
D   115.02694
E   129.04259
F   147.06841
G   57.02146
H   137.05891
I   113.08406
K   128.09496
L   113.08406
M   131.04049
N   114.04293
P   97.05276
Q   128.05858
R   156.10111
S   87.03203
T   101.04768
V   99.06841
W   186.07931
Y   163.06333"
    |> fun x->x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)
    |> Seq.map (fun x-> x.Split([|' '|],StringSplitOptions.RemoveEmptyEntries))
    |> Seq.map (fun toks-> (toks.[0],Double.Parse(toks.[1])))
    |> Map.ofSeq

let monomass str =
    str
    |> Seq.map (fun c-> monomasstbl.[string c])
    |> Seq.sum

printfn "%f" (monomass "SKADYEK")      

printfn "%f" (monomass "FRKYGVCLEKERMHQAKFGNLQYPLSMLTWKHDIPTTMDRTMIIKQAATDMCDYGNTSYVHVGASWDIMNRETKHVYVFFSVTQVINDNNDAHWNDVGHVLEWPAKRYCQHTQHWCEAHCMPEAMSHHAITWHTDHNHNFTHHCGQGTSFICEFAQGNNLEGHVDFHRIAELARTIMPGYNDCENDKEHCPPWPWFINRTLMDTWWVNVIREFKYEAIRMLFPQKMYWHMNCRCHPKGTYRIHKLPQHGHMSYKRMGDMGIEILTIEEPWTNGIWCVQDEPDEIIIMMGCMAPPYTQQDPECFEDAKWCFWQRQRNTVHVGRMTKVFNNMVMIHQIRHCFHFIPDGLGDAPLAVYVNETEPANVKGGSKTMPVPPFTHSRCSALAHRTLQNTHTAEWASWDARARFVCEWHRVPNRFCGHNYGQMGPLHHRVFNNTEMHRVTDQKVVNYPETSQTFKFCWGETHCMNFGQCLSPEGVICHFALPSYEICSGRFMDYKVGAVVPDMHFAAVYKWKDCGFWNCQMRYRPRTKYVAEDWFAVPWRMASVVLPHRINTIVCWALGDHEVQNQLDRCIFRGTNQCAPIAGDQWFGNHKHTWPKSHTVNVQMVGIGEGVESILYLSTILSEPFGFIQSMMVPHPVVVPWTCAMNCLAYEYQIVEGIPWKPCGCYNCIWPTKNGLQYAGAATQTMVRDEMAQFCFTYPHWIKGFYGVKKQNPHSTRMVLRGFNMSQRLEGCTVAYYPCESEEFVHHYHAFNCPMWDWWMTTWHDDYKMLCLTHAHRCCFRTTTFQWAVYILQWQDQMNANQLSPEKMWSPTAPVYHQLREVLRMGGCFYAYTPFPGPHIDHLAYVETQCYH")

// revp
let isRevPalindrome (str:string) = // probably not very efficient.
    str=reverseComplement str

let enumerateSomeSubstrings min max (str:string) =
    seq {
        for i in [0..str.Length-min] do
            for j in [min..max] do
                if (i+j<=str.Length) then
                    yield (i,j,str.Substring(i,j))
    }

getData "revp"
|> fun x->x.Trim()
|> parseFasta
|> Seq.exactlyOne
|> fun x->enumerateSomeSubstrings 4 12 x.String
//|> Seq.iter (printfn "%A")
|> Seq.filter (fun (i,j,s)->isRevPalindrome s)
|> Seq.iter (fun (i,j,_) -> printfn "%d %d" (i+1) j)   


// splc
let splice (dnastr:FASTA) (introns:FASTA list) =
    introns
    |> Seq.fold (fun (s:string) (i:FASTA) -> s.Replace(i.String,"")) dnastr.String
    
// According to rosalind there should be only one solution but parseProteins returns several.  I took the longest.
getData "splc"
|> fun x -> x.Trim()
|> parseFasta
|> List.ofSeq
|> fun l -> 
    match l with 
    | head::tail -> splice head tail
    | _ -> ""
|> dnaToRna
|> parseProteins
|> Seq.iter (printfn "%s")

// lexf
getData "lexf"
|> splitNewline
|> fun toks -> (toks.[0].Split(' ')|>Array.toSeq,Int32.Parse(toks.[1]))
|> fun (alphabet:seq<string>,strlen:int) ->
    let rec addChar (i:int) (s:string) : seq<string> =
        if i=0 then
            seq {yield s}
        else
            seq {               
                for a in alphabet do
                    yield! (s+a)|>addChar (i-1)
            }
    addChar strlen ""
|> Seq.iter (printfn "%s")

// lgis. unfinished.
// This was the first attempt. Too expensive.
let longestSubsequence2 cmp s =
    Seq.fold (fun (sequences:list<list<int>>) (elem:int) ->
                sequences
                |> List.map 
                    (fun si ->
                        match si with
                        | head :: tail -> 
                            if cmp elem head then
                                [elem::si;si]
                            else
                                [si]
                        | [] -> []
                    )
                |> List.concat
                |> fun l -> [elem]::l
                ) 
                [] s
    |> List.maxBy (fun s -> List.length s)
    |> List.rev // needs reversing when done.

// On the second attempt we build a data structure: what is the longest subsequence that ends with p<x?
// this is a question of building a sorted list: a red-black tree will do fine and stackoverflow has a nice implementation here:
// http://stackoverflow.com/questions/20297431/difficulty-in-writing-red-black-tree-in-f
// This will work but it's likely still too expensive...
type color = R | B
type 'a tree = E | T of color * 'a tree * 'a * 'a tree

let balance = function
  | B, T (R, T (R,a,x,b), y, c), z, d
  | B, T (R, a, x, T (R,b,y,c)), z, d
  | B, a, x, T (R, T (R,b,y,c), z, d)
  | B, a, x, T (R, b, y, T (R,c,z,d)) -> T (R, T (B,a,x,b), y, T (B,c,z,d))
  | col, a, x, b                      -> T (col, a, x, b) 

let insert cmp x s = 
  let rec ins = function
    | E                  -> T (R,E,x,E)
    | T (col,a,y,b) as s ->
        if x=y then
          s
        elif cmp x y then
          balance (col, ins a, y, b)
        else 
          balance (col, a, y, ins b)
        
  match ins s with
  | T (_,a,y,b) -> T (B,a,y,b)
  | t -> t

// This is an inorder traversal of all nodes such that x cmp node
let rec leftLookup t cmp x : 'a list =    
    match t with
    | E -> []
    | T(_,left,curr,right) ->
        match x with // when passed None we return in-order traversal of the whole tree.
        | None -> [leftLookup left cmp x;[curr];leftLookup right cmp x] |> List.concat
        | Some(x') -> // otherwi
            if cmp x' curr then
                leftLookup left cmp x
            else 
                [leftLookup left cmp x;[curr];leftLookup right cmp x]
                |> List.concat

// unroll the prevls. see below. in progress.
let rec lookup t cmp x =
    match t with
    | E -> None
    | T (_,left,(curr,prev,prevl),right) ->
        if x=curr then Some (curr,prev,prevl)
        else if cmp x curr then lookup left cmp x
        else lookup right cmp x

let unroll (t:(int*int*int) tree) cmp lastElt =
    Seq.unfold 
        (fun prev -> 
            let pnext=lookup t cmp prev
            match pnext with
            | None -> None
            | Some (_,pp,_) -> Some (prev,pp)
        ) lastElt

let longestSubsequence cmp s =
    let cmp' (x,_,_) (y,_,_) = cmp x y
    Seq.fold (fun (t: (int*int*int) tree) (elem:int) ->
        leftLookup t cmp' (Some (elem,-1,-1)) // find all sequences that end with values < elem
        |> fun l -> 
            if List.isEmpty l then 
                insert cmp' (elem,-1,1) t
            else
                l
                |> Seq.maxBy (fun (x,prev,prevl) -> prevl)
                |> fun (x,prev,prevl) -> insert cmp' (elem,x,prevl+1) t) E s
    |> fun t -> // now we have the tree.  how do we get the largest prevl?
        leftLookup t cmp' None
        |> List.maxBy (fun (x,prev,len) -> len)
        |> fun (lastElt,_,_) -> unroll t cmp lastElt
    |> List.ofSeq
    |> List.rev    
         


getData "lgis"
|> splitNewline
|> fun toks -> toks.[1]
|> fun x -> x.Split(' ')
|> Array.map (fun x -> Int32.Parse(x))
|> fun x -> x,longestSubsequence (<) x
|> fun (x,l1) -> (x,l1,longestSubsequence (>) x)
|> fun (x,l1,l2) ->
    printfn "%s" (l1|>List.map string|> String.concat " ")
    printfn "%s" (l2|>List.map string|> String.concat " ")

// TRAN
type Mutation = TRANSITION | TRANSVERSION | NONE
let whatMutation (c1,c2) =
    match c1 with
    | 'A' -> 
        if c2='G' then TRANSITION
        else if c2='A' then NONE
        else TRANSVERSION
    | 'C' ->
        if c2='T' then TRANSITION
        else if c2='C' then NONE
        else TRANSVERSION
    | 'G' ->
        if c2='A' then TRANSITION
        else if c2='G' then NONE
        else TRANSVERSION
    | 'T' ->
        if c2='C' then TRANSITION
        else if c2='T' then NONE
        else TRANSVERSION
    | _ -> raise (new Exception("Unrecognized nucleotide"))

let ttRatio s1 s2 =
    Seq.zip s1 s2
    |> Seq.map whatMutation  
    |> Seq.countBy id  
    |> Map.ofSeq
    |> fun m -> ((float m.[TRANSITION])/(float m.[TRANSVERSION]))

getData "tran"
|> fun x->x.Trim()
|> parseFasta
|> Array.ofSeq
|> fun t -> ttRatio t.[0].String t.[1].String
|> printfn "%f"


// LONG
type SuperstringMatch =
    Superstring of string | Disjoint of (string*string) 
// Some things to improve: return option instead of disjoint.
let combineStrings (s1:string) (s2:string) =        
    let (a,b)=if s1.Length<s2.Length then (s1,s2) else (s2,s1)
    [-a.Length .. b.Length ]
    |> Seq.map 
        (fun i -> 
        //    printfn "i: %d" i
            if i<0 then
                let overlap= a.Length+i
      //          printfn "o: %d" overlap
                if overlap=0 then Disjoint(a,b)
                else                                        
                    let b'=b.Substring(0,overlap)
                    let a'=a.Substring(-i,overlap)
                    if a'=b' then Superstring(a+(b.Substring(overlap)))
                    else Disjoint (a,b)
            else
                if (i<b.Length-a.Length) then
                    let overlap=a.Length
                    let b'=b.Substring(i,overlap)
                    if b'=a then Superstring(b)
                    else Disjoint(a,b) // obviously wrong
                else 
                    let overlap = b.Length-i
    //                printfn "o2: %d" overlap
                    if overlap=0 then Disjoint(a,b)
                    else
                        let b'=b.Substring(i,overlap)
                        let a'=a.Substring(0,overlap)
                        if a'=b' then Superstring(s1+(s2.Substring(overlap)))
                        else Disjoint(a,b))
    |> Seq.countBy id
    |> Seq.map (fun (a,b)->a)

let cost ssm = 
    match ssm with
    | Disjoint(a,b) -> a.Length+b.Length
    | Superstring(s) -> s.Length


let rec mergeSS s1 s2 =
    match s1 with
    | Disjoint (a,b) ->
        match s2 with 
        | Disjoint (c,d) ->
            [s1;s2]
        | Superstring(ss2) -> 
            [combineStrings a ss2; combineStrings b ss2]|> Seq.concat|>List.ofSeq           
    | Superstring(ss1) ->
        match s2 with
        | Disjoint (c,d) ->
            [combineStrings c ss1; combineStrings d ss1]|> Seq.concat|>List.ofSeq           
        | Superstring (ss2) ->
            combineStrings ss1 ss2 |> List.ofSeq    


let findBest ssm s =
    s
    |> Seq.map (fun x -> (x,mergeSS ssm (Superstring(x))|>Seq.minBy cost))        
    |> Seq.minBy (fun (x,s) -> cost s)

// to improve perf: switch to set indexing on labels rather than full identity.  Needs an extra map.
let findSuperstring l : SuperstringMatch =
    match l with 
    | head::tail -> 
        let (startx,startset)=Superstring(head),tail|>Set.ofList
        
        Seq.unfold 
            (fun (currstring,currset) -> 
                if (Set.isEmpty currset) then 
                    None
                else
                    let (remove,nextString)=findBest currstring currset
                    Some(nextString,(nextString,Set.remove remove currset))                    
                    ) (startx,startset)
        |> List.ofSeq
        |> List.rev
        |> Seq.take 1
        |> Seq.exactlyOne
    | [] -> Superstring("")

getData "long"
|> fun x -> x.Trim()
|> parseFasta
|> Seq.map (fun x -> x.String)
|> List.ofSeq
|> findSuperstring
|> printfn "%A"

// PMCH
// pn is basically n!  Implemented with bigint since we overlow int64.
let rec pn (n:bigint) =
    if n=bigint 1 then bigint 1
    else n*(pn (n-bigint 1)) // only n ways to connect since the graph must be bipartite. This is a slightly different scenario from what's described in paragraph 2 of the problem description.

let pnAUCG (auCount:int,cgCount:int) =
    (pn (bigint auCount))*(pn (bigint cgCount))
//    (pn auCount)*(pn cgCount)

getData "pmch"
|> fun x -> x.Trim()
|> parseFasta
|> Seq.exactlyOne
|> fun x-> x.String|>Seq.countBy id
|> Map.ofSeq
|> fun m -> (m.['A'],m.['C']) // counting pairs so only take one side of each pair in the count
|> pnAUCG
|> printfn "%A"

// PPER
let rec P n k =
    if k=1 then n
    else (n*(P (n-1) (k-1)))%1000000

// input: 
(*
> P 21 7;;
val it : int = 51200
> *)

// PROB
let pRand inStr gcContent =
    [('A',(1.-gcContent)/2.);('C',gcContent/2.);('G',gcContent/2.);('T',(1.-gcContent)/2.)]
    |> Map.ofList
    |> fun m -> (inStr |> Seq.map (fun c -> m.[c]|>log10) |> Seq.sum)

getData "prob"
|> splitNewline
|> fun toks -> toks.[0], toks.[1].Split(' ')|>Seq.map (Double.Parse)
|> fun (inStr,gcContents) -> (gcContents|>Seq.map (pRand inStr))
|> Seq.map (sprintf "%f")
|> String.concat " "
|> printfn "%s"


// SIGN
let enumerate N = 
    [0..int (2.**(float N)-1.)]
    |> Seq.map (
        fun X ->
            Seq.unfold 
                (fun (x,i) -> 
                    if i=0 then None 
                    else
                        if x%2=0 then Some (1,(x/2,i-1)) 
                        else Some (-1,(x/2,i-1))) (X,N)|>Array.ofSeq)

let cross f s1 s2 =
    seq {
        for a in s1 do
            yield! seq {
                for b in s2 do
                    yield f a b
                }
        }

let applySign p s =
    Seq.zip p s
    |> Seq.map (fun (a,b) -> a*b)

let printVector v =
    v
    |> Seq.map string
    |> String.concat " "
    |> printfn "%s"

let doSign N = 
    let result = cross applySign (enumerate N) (permute [1..N]) |> List.ofSeq
    printfn "%d" (List.length result)
    result
    |>Seq.iter printVector
    

// SSEQ
let findSseq (s:string) (t:string) =
    Seq.unfold (fun (i,j) ->
        if (i=s.Length || j=t.Length) then None 
        else
            let ix = s.IndexOf(t.[j],i)
            Some (ix,(ix+1,j+1))) (0,0)
    |> Seq.map ((+) 1)

getData "sseq"
|> fun x -> x.Trim()
|> parseFasta
|> Array.ofSeq
|> fun toks -> findSseq (toks.[0].String) (toks.[1].String)
|> Seq.map string
|> String.concat " "
|> printfn "%s"

// TREE
// we'll assume all the inputs are trees.  We just need to count CCs.
let countCCs vCount edges =
    let initLookup = [1..vCount]|> Seq.map (fun x ->x,x)|> Map.ofSeq
    let initSets = [1..vCount]|> Seq.map (fun x -> (x,[x]))|>Map.ofSeq
    edges
    |>Seq.fold 
        (fun (componentLookup:Map<int,int>,sets:Map<int,int list>,nextIx) (v1,v2) -> 
            if Map.containsKey v1 componentLookup then
                let c1=componentLookup.[v1]
                if Map.containsKey v2 componentLookup then
                    let c2=componentLookup.[v2]
                    if c1=c2 then (componentLookup,sets,nextIx)
                    else
                        if List.length sets.[c2]>List.length sets.[c1] then
                            let outLookup = Seq.fold (fun m vi -> Map.add vi c2 m) componentLookup sets.[c1]
                            let outSets = sets|>Map.add c2 (List.concat [sets.[c2];sets.[c1]]) |> Map.remove c1
                            (outLookup,outSets,nextIx)
                        else 
                            let outLookup = Seq.fold (fun m vi -> Map.add vi c1 m) componentLookup sets.[c2]
                            let outSets = sets|>Map.add c1 (List.concat [sets.[c1];sets.[c2]]) |> Map.remove c2
                            (outLookup,outSets,nextIx)
                else
                    (Map.add v2 c1 componentLookup, Map.add c1 (v2::sets.[c1]) sets, nextIx)
            else
                if Map.containsKey v2 componentLookup then
                    let c2=componentLookup.[v2]
                    (Map.add v1 c2 componentLookup, Map.add c2 (v1::sets.[c2]) sets, nextIx)
                else
                    let newSet=nextIx
                    (Map.add v1 newSet componentLookup|>Map.add v2 newSet, Map.add newSet [v1;v2] sets, newSet+1)
    ) (initLookup,initSets,vCount+1)
    |> (fun (_,s,_) -> s |> Map.toList |> List.length) 
   
getData "tree"
|> fun x -> x.Trim()
|> splitNewline
|> fun toks ->
    (int toks.[0], toks |>Seq.skip 1 |> Seq.map (fun s -> s.Split(' ') |> fun t -> int t.[0],int t.[1]))
|> fun (count, edges) -> countCCs count edges
|> fun x -> x-1   
|> printfn "%d"

// CAT
// first let's tool with Catalan numbers.  They won't quite represent the answer because we have additional constraints due to base pair names.
// this is an expensive approach since it fails to capture partial results.
// How to do proper dynamic programming in a functional way?  With a mutable map?

// In the end the solution is very easy: there are a limited number of split points which have equal 
// pairings on both the left and the right.  Find those splits and count up. function could only be 
// improved with proper dynamic accounting.
open System.Collections.Generic
 // some handy dynamic programming snippets from http://www.fssnip.net/8P
 /// The function creates a function that calls the argument 'f'
 /// only once and stores the result in a mutable dictionary (cache)
 /// Repeated calls to the resulting function return cached values.
let memoize f =    
   // Create (mutable) cache that is used for storing results of 
   // for function arguments that were already calculated.
    let cache = new Dictionary<_, _>()
    (fun x ->
       // The returned function first performs a cache lookup
       let succ, v = cache.TryGetValue(x)
       if succ then v else 
         // If value was not found, calculate & cache it
         let v = f(x) 
         cache.Add(x, v)
         v)

let isComplement c1 c2 =
    (c1='A' && c2='U')
    || (c1='U' && c2='A')
    || (c1='C' && c2='G')
    || (c1='G' && c2='C')

let splitStr (str:string) (x:int) =
    // split the string at x, dropping the first char, and the char at x
    (str.Substring(1,x-1),str.Substring(x+1,str.Length-x-1))

let isValidStr = memoize (fun str ->
    let m=
        str
        |> Seq.countBy id
        |> Map.ofSeq
    (guard m 'A')=(guard m 'U') && (guard m 'C')=(guard m 'G'))

let isValidSplit (l,r) =
    isValidStr l && isValidStr r

let rec cat = memoize (fun str ->     
    if not (isValidStr str) then raise (new Exception (sprintf "Invalid string %s" str))
    if str="" || str.Length=2 then 1L    
    else
        [1..str.Length-1] // this can be sped up by jumping in increments of 2.
        |> Seq.map 
            (fun x -> 
                if isComplement str.[0] str.[x] then 
                    let (l,r)=splitStr str x
                    if isValidSplit (l,r) then
                        (cat l * cat r) 
                    else
                        0L
                else 0L)
        |> Seq.sum
        |> fun x -> x%1000000L)
    
getData "cat"
|> parseFasta
|> Seq.exactlyOne
|> fun x -> cat x.String
|> printfn "%d"    
    
// corr:
// We'll go with the constraints rather than consensus
let corrData=
    getData "corr"
    |> parseFasta
    |> Seq.map (fun x -> [x.String;reverseComplement x.String])
    |> List.concat
    |> Seq.countBy id
    |> Map.ofSeq

let correctReads=
    corrData
    |> Map.toSeq
    |> Seq.map (fun (x,c) -> (x,c+(guard corrData (reverseComplement x))))
    |> Seq.filter (fun (x,c)->c>=4) // if it appears twice in the data set we'll have also replicated it twice.
    |> List.ofSeq

let closestRead x =
    correctReads
    |> Seq.map (fun (y,c)-> (y,dh x y ))
    |> Seq.minBy (fun (y,hamm)-> hamm)
    

getData "corr"
|> parseFasta
|> Seq.map (fun x -> (x.String,closestRead x.String))
|> Seq.filter (fun (x,(y,h))-> h=1)
|> Seq.map (fun (x,(y,h)) -> sprintf "%s->%s" x y)
|> Seq.iter (printfn "%s")    



// inod
// wikipedia helps here. https://en.wikipedia.org/wiki/Unrooted_binary_tree
// An unrooted tree with 3 leaves has 1 internal node.
// adding another leaf requires splitting and edge and adding a vertex.  So a tree with 4 leaves has 2 internal nodes.
// in general n_l = l-2
// 
getData "inod"
|> int
|> fun x -> x-2
|> printfn "%d"

// kmer. Would be trivial if we didn't have to enumerate all possible 
let rec enumerateKmers k : string list =
    if k=0 then [""]
    else 
        enumerateKmers (k-1)
        |> List.map (fun l ->
            [
                "A"+ l 
                "C" + l
                "G" + l
                "T" + l
            ])
        |> List.concat
            

      
getData "kmer"
|> parseFasta
|> Seq.exactlyOne
|> fun x -> x.String
|> Seq.windowed 4
|> Seq.map (fun c -> c|>Seq.map string|>String.concat "")
|> Seq.fold (fun m s-> Map.add (string s) (m.[string s]+1) m) (enumerateKmers 4|>Seq.map (fun s -> (s,0))|> Map.ofSeq)
|> Map.toSeq
|> Seq.sortBy (fun (x,y)->x)
|> Seq.map (fun (x,y)-> string y)
|> String.concat " "
|> printfn "%s"


// KMP

// Failure array for Knuth-Morris-Pratt
// https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm
// a bit of a hack adding List<int>, so we have O(1) lookback.
// This did not work.  Somewhere there is an off-by-one issue.
//open System.Collections.Generic
(* not quite right.
let failureArray (s:string) =
    let mutable t=Array.zeroCreate s.Length
    Seq.unfold (
        fun (pos:int,cnd:int) ->
            eprintfn "%d\t%d" pos cnd
            if pos=s.Length then None
            else if s.[pos] = s.[cnd] then                    
                t.[pos] <- cnd + 1
                Some (cnd+1,(pos+1, (cnd+1)))
            else if cnd > 0 then                
                Some (t.[cnd],(pos, t.[cnd]))
            else 
                t.[pos] <- 0 
                //if cnd>0 then eprintfn "%s" "unexpected cnd"
                Some (cnd,(pos+1,cnd)))
        (2,0)
    |> fun s -> (s|>List.ofSeq,t)    
    *)
// http://www.inf.fh-flensburg.de/lang/algorithmen/pattern/kmpen.htm
// differs slightly from the web page by offsetting the boundary array by 1.
// This is probably the one problem in this set where laziness got the best of me and I went with the online code.
// I fundamentally detest array twiddling- glad other smart people work these things out.
let failureArray (p:string) =
    let mutable i=0
    let mutable j= -1
    let mutable b=Array.zeroCreate (p.Length+1)
    b.[i]<-j;
    while (i<p.Length) do
        eprintfn "%d\t%d" i j
        while (j>=0 && p.[i]<>p.[j]) do
            j<-b.[j]
        i<-i+1 
        j<-j+1
        b.[i]<-j;
    b
    |> Seq.ofArray
    |> Seq.skip 1
    

getData "kmp"
|> parseFasta
|> Seq.exactlyOne
|> fun x-> failureArray x.String
|> Seq.map string
//|> fun (s,t) -> t|>Array.map string
|> String.concat " "
|> printfn "%s"

// LCSQ
// Greedy approach.
let lcsq (s1:string) (s2:string) =
    let rec lcsqhelper = 
        memoize (fun (i,j) ->
        if i=s1.Length || j=s2.Length then ""
        else
            if s1.[i]=s2.[j] then 
                string s1.[i]+(lcsqhelper (i+1,j+1))            
            else
                // should be possible to prune certain branches but I ended up with the wrong answer when I tried this. Dynamic programming FTW.
                let sq1=lcsqhelper(i+1,j)
                let sq2=lcsqhelper(i,j+1)
                if sq1.Length>sq2.Length then
                    sq1
                else
                    sq2
     )
    lcsqhelper (0,0)       

let rec isSubsequence (s:string) (sq:string) =
    if sq.Length=0 then true
    else if s.Length=0 then false
    else if s.[0]=sq.[0] then isSubsequence (s.Substring(1)) (sq.Substring(1))
    else isSubsequence (s.Substring(1)) sq
    

getData "lcsq"
|> parseFasta
|> Array.ofSeq
|> fun toks -> (toks.[0].String,toks.[1].String,lcsq toks.[0].String toks.[1].String)
|> fun (s1,s2,sq) -> (isSubsequence s1 sq),(isSubsequence s2 sq),sq
|> fun (b1,b2,sq) -> printfn "%b\n%b\n%s\n%d" b1 b2 sq sq.Length

// lexv
// variation on lexf.
getData "lexv"
|> splitNewline
|> fun toks -> (toks.[0].Split(' ')|>Array.toSeq,Int32.Parse(toks.[1]))
|> fun (alphabet:seq<string>,strlen:int) ->
    let alphabet' = "_" :: (List.ofSeq alphabet)
    let rec addChar (i:int) (s:string) : seq<string> =
        if i=0 then
            seq {yield s}
        else
            seq {               
                for a in alphabet' do
                    if a="_" then yield s
                    else 
                        yield! (s+a)|>addChar (i-1)
            }
    addChar strlen ""
    //|> Seq.map (fun s->s.Replace("_",""))
    |> Seq.filter (fun s-> not (String.IsNullOrEmpty(s)))
|> Seq.iter (printfn "%s")


[<EntryPoint>]
let main argv = 
    // dna
    File.ReadAllText(argv.[0])
    |> fun x -> x.Trim()
    |> Seq.countBy (fun x->x)
    |> Map.ofSeq
    |> fun m -> printfn "%d %d %d %d" (guard m 'A') (guard m 'C') (guard m 'G') (guard m 'T')
    0 // return an integer exit code




