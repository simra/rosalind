// Learn more about F# at http://fsharp.net
// See the 'F# Tutorial' project for more help.
open System
open System.IO

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)

// 1. DNA
let guard m c =
    if Map.containsKey c m then
        m.[c]
    else
        0

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
type FASTA =
    {Label:string;String:string}
    
let parseFasta (text:string) =
    text.Split([|'>'|],StringSplitOptions.RemoveEmptyEntries)
    |> Seq.map 
        (fun lines-> 
            lines.Split([|'\n';'\r'|],StringSplitOptions.RemoveEmptyEntries)
            |> fun toks -> 
                {Label=toks.[0]; String=toks.[1..]|>String.concat ""})

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

let signPerm n =
    perm n



[<EntryPoint>]
let main argv = 
    // dna
    File.ReadAllText(argv.[0])
    |> fun x -> x.Trim()
    |> Seq.countBy (fun x->x)
    |> Map.ofSeq
    |> fun m -> printfn "%d %d %d %d" (guard m 'A') (guard m 'C') (guard m 'G') (guard m 'T')
    0 // return an integer exit code




