/** -*-c++-*-
 *
 *  Copyright 2014 The Australian National University
 *
 *  ImageAnalysis.hpp
 *
 *  Facade classes for volume image data and associated information such as
 *  Morse vector field, Morse chain complex and cell labellings.
 *
 *  Olaf Delgado-Friedrichs mar 14
 *
 */

#ifndef ANU_AM_DIAMORSE_IMAGEANALYSIS_HPP
#define ANU_AM_DIAMORSE_IMAGEANALYSIS_HPP


#include <cmath>
#include <string>
#include <vector>

#include "callables.hpp"
#include "chainComplexExtraction.hpp"
#include "CubicalComplex.hpp"
#include "MorseVectorField.hpp"
#include "PackedMap.hpp"
#include "Partition.hpp"
#include "persistence.hpp"
#include "SimpleComplex.hpp"
#include "simplification.hpp"
#include "traversals.hpp"
#include "Set.hpp"
#include "vectorFieldExtraction.hpp"
#include "VertexMap.hpp"
#include "volume_io.hpp"


namespace anu_am
{
namespace diamorse
{


template<typename Value> class MorseData;


template<typename Value>
class ImageData
{
    typedef CubicalComplex::cell_id_type Cell;
    typedef VertexMap<CubicalComplex, Value> Scalars;

    std::vector<size_t> dims_;
    CubicalComplex complex_;
    Scalars scalars_;

    friend class MorseData<Value>;

public:
    ImageData();

    ImageData(std::string const filename)
        : dims_(readDimensions(filename)),
          complex_(dims_.at(0), dims_.at(1), dims_.at(2)),
          scalars_(complex_, readVolumeData<Value>(filename))
    {
    }

    int xdim() const
    {
        return complex_.xdim();
    }

    int ydim() const
    {
        return complex_.ydim();
    }

    int zdim() const
    {
        return complex_.zdim();
    }

    std::vector<float> positionForCell(Cell const v) const
    {
        return complex_.cellPosition(v);
    }

    Cell cellAtPosition(std::vector<float> const p) const
    {
        return complex_.cellAt(p.at(0), p.at(1), p.at(2));
    }

    int cellDimension(Cell const v) const
    {
        return complex_.cellDimension(v);
    }

    Value cellValue(Cell const v) const
    {
        Vertices vertices(complex_.xdim(), complex_.ydim(), complex_.zdim());

        Value val = scalars_.get(vertices(v, 0));
        for (size_t i = 1; i < vertices.count(v); ++i)
            val = std::max(val, scalars_.get(vertices(v, i)));

        return val;
    }
};


template<typename Value>
class MorseData
{
    typedef CubicalComplex::cell_id_type Cell;
    typedef MorseVectorField<PackedMap> Field;
    typedef Field::DataPtr::element_type::value_type FieldItem;
    typedef VertexMap<CubicalComplex, Cell> Labels;
    typedef VertexMap<CubicalComplex, float> Depths;
    typedef VertexMap<CubicalComplex, bool> Mask;
    typedef std::vector<std::pair<Cell, int> > Boundary;
    typedef std::map<Cell, Boundary> MorseComplex;

    struct Comparator
    {
        Comparator(ImageData<Value> const& img)
            : img_(img)
        {
        }

        bool operator()(Cell const v, Cell const w)
        {
            return img_.cellValue(v) < img_.cellValue(w)
                or (img_.cellValue(v) == img_.cellValue(w)
                    and img_.cellDimension(v) < img_.cellDimension(w));
        }

    private:
        ImageData<Value> const& img_;
    };

    ImageData<Value> imageData_;
    float threshold_;

    mutable Field field_;
    mutable bool fieldComputed_;

    mutable std::vector<Cell> critical_;
    mutable bool criticalComputed_;

    mutable std::vector<Cell> criticalSorted_;
    mutable bool criticalSortedComputed_;

    mutable Labels basins_;
    mutable bool basinsComputed_;

    mutable Mask watersheds_;
    mutable bool watershedsComputed_;

    mutable Depths skeleton_;
    mutable bool skeletonComputed_;

    mutable Mask paths_;
    mutable bool pathsComputed_;

    mutable MorseComplex chains_;
    mutable bool chainsComputed_;

    mutable MorseComplex cochains_;
    mutable bool cochainsComputed_;

    mutable SimpleComplex simpleChainComplex_;
    mutable bool simpleChainComplexComputed_;

    mutable std::vector<std::pair<Cell, Cell> > birthDeathPairs_;
    mutable bool birthDeathPairsComputed_;

    void runSimplify() const
    {
        VertexMap<CubicalComplex, Value> const scalars = imageData_.scalars_;
        CubicalComplex complex = imageData_.complex_;
        Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());

        simplify(complex,
                 field_,
                 withArgument(maxima(scalars, vertices)),
                 mayCancel(complex, scalars, field_, threshold_));
    }

    Field const& field() const
    {
        if (not fieldComputed_)
        { 
            fillMorseVectorField(
                imageData_.complex_,
                imageData_.scalars_,
                field_);

            if (threshold_ >= 0)
                runSimplify();

            fieldComputed_ = true;
        }
        return field_;
    }

    std::vector<Cell> const& critical() const
    {
        if (not criticalComputed_)
        {
            CubicalComplex complex = imageData_.complex_;

            for (Cell cell = 0; cell <= complex.cellIdLimit(); ++cell)
                if (complex.isCell(cell) and field().isCritical(cell))
                    critical_.push_back(cell);

            criticalComputed_ = true;
        }
        return critical_;
    }

    std::vector<Cell> const& criticalSorted() const
    {
        if (not criticalSortedComputed_)
        {
            criticalSorted_ = critical();
            std::stable_sort(criticalSorted_.begin(), criticalSorted_.end(),
                             Comparator(imageData_));

            criticalSortedComputed_ = true;
        }
        return criticalSorted_;
    }

    Labels const& basins() const
    {
        if (not basinsComputed_)
        {
            basins_ = Labels(imageData_.complex_,
                             imageData_.complex_.cellIdLimit() + 1);
            CubicalComplex const& complex = imageData_.complex_;
            Field::Vectors coV = field().coV();
            Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);

            std::vector<Cell> const& sources = critical();

            for (size_t i = 0; i < sources.size(); ++i)
            {
                Cell const s = sources.at(i);
                if (complex.cellDimension(s) != 0)
                    continue;
                basins_.set(s, s);

                std::vector<std::pair<Cell, Cell> > t =
                    flowTraversal(s, coV, coI);
                for (size_t j = 0; j < t.size(); ++j)
                {
                    Cell const b = t.at(j).second;
                    Cell const c = coV(b);
                    if (b != c)
                        basins_.set(c, s);
                }
            }

            basinsComputed_ = true;
        }

        return basins_;
    }

    Mask const& watersheds() const
    {
        if (not watershedsComputed_)
        {
            watersheds_ = Mask(imageData_.complex_, false);
            CubicalComplex const& complex = imageData_.complex_;
            Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());
            Field::Vectors coV = field().coV();
            Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);

            std::vector<Cell> const& sources = critical();

            for (size_t i = 0; i < sources.size(); ++i)
            {
                Cell const s = sources.at(i);
                if (complex.cellDimension(s) != 1)
                    continue;

                std::vector<std::pair<Cell, Cell> > t =
                    flowTraversal(s, coV, coI);
                for (size_t j = 0; j < t.size(); ++j)
                {
                    Cell const b = t.at(j).second;
                    Cell const c = coV(b);
                    if (b != c) {
                        size_t const m = vertices.count(c);
                        for (size_t i = 0; i < m; ++i)
                            watersheds_.set(vertices(c, i), true);
                    }
                }
            }

            watershedsComputed_ = true;
        }
        return watersheds_;
    }

    Depths const& skeleton() const
    {
        if (not skeletonComputed_)
        {
            skeleton_ = Depths(imageData_.complex_, INFINITY);
            Field::Vectors V = field().V();
            CubicalComplex const& complex = imageData_.complex_;
            Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
            Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());

            std::vector<Cell> const& sources = critical();

            for (size_t k = 0; k < sources.size(); ++k)
            {
                Cell const s = sources.at(k);
                Value const val = imageData_.cellValue(s);

                std::vector<std::pair<Cell, Cell> > t =
                    flowTraversal(s, V, I);
                for (size_t j = 0; j < t.size(); ++j)
                {
                    Cell const c = V(t.at(j).second);

                    size_t const m = vertices.count(c);
                    for (size_t i = 0; i < m; ++i)
                    {
                        Cell const v = vertices(c, i);
                        skeleton_.set(v, std::min(skeleton_.get(v), val));
                    }
                }
            }

            skeletonComputed_ = true;
        }

        return skeleton_;
    }

    Mask const& paths() const
    {
        if (not pathsComputed_)
        {
            paths_ = Mask(imageData_.complex_, false);
            Field::Vectors V = field().V();
            Field::Vectors coV = field().coV();
            CubicalComplex const& complex = imageData_.complex_;
            Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
            Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);
            Vertices vertices(complex.xdim(), complex.ydim(), complex.zdim());

            std::vector<Cell> const& sources = critical();

            std::vector<bool> active(imageData_.complex_.cellIdLimit() + 1);
            std::queue<Cell> queue;

            for (size_t i = 0; i < sources.size(); ++i)
                queue.push(sources.at(i));

            while (not queue.empty())
            {
                Cell const a = queue.front();
                queue.pop();

                size_t const n = coI.count(a);
                for (size_t i = 0; i < n; ++i)
                {
                    Cell const b = coI(a, i);

                    if (coV.defined(b) and coV(b) != a )
                    {
                        Cell const c = coV(b);

                        if (c != b and not active.at(c))
                        {
                            active.at(c) = true;
                            queue.push(c);
                        }
                    }
                }
            }

            std::set<Cell> seen;
            for (size_t i = 0; i < sources.size(); ++i)
                queue.push(sources.at(i));

            while (not queue.empty())
            {
                Cell const a = queue.front();
                queue.pop();

                size_t const n = I.count(a);
                for (size_t i = 0; i < n; ++i)
                {
                    Cell const b = I(a, i);

                    if (V.defined(b) and V(b) != a )
                    {
                        Cell const c = V(b);

                        if (c != b and active.at(b) and seen.count(c) == 0)
                        {
                            seen.insert(c);
                            for (size_t j = 0; j < vertices.count(c); ++j)
                                paths_.set(vertices(c, j), true);
                            queue.push(c);
                        }
                    }
                }
            }

            pathsComputed_ = true;
        }
        return paths_;
    }

    MorseComplex const& chains() const
    {
        if (not chainsComputed_)
        {
            chains_ = chainComplex(imageData_.complex_, field());
            chainsComputed_ = true;
        }
        return chains_;
    }

    MorseComplex const& cochains() const
    {
        if (not cochainsComputed_)
        {
            cochains_ = reverse(chains());
            cochainsComputed_ = true;
        }
        return cochains_;
    }

    SimpleComplex const& simpleChainComplex() const
    {
        if (not simpleChainComplexComputed_)
        {
            std::vector<Cell> const& sources = criticalSorted();
            size_t const n = sources.size();

            std::map<Cell, size_t> index;
            for (size_t i = 0; i < n; ++i)
                index[sources.at(i)] = i;

            std::vector<unsigned int> dims;
            std::vector<float> scalars;
            std::vector<std::vector<Cell> > faceLists;

            for (size_t i = 0; i < n; ++i)
            {
                Cell const v = sources.at(i);
                dims.push_back(imageData_.cellDimension(v));
                scalars.push_back(imageData_.cellValue(v));

                std::vector<Cell> const flin = chainComplexBoundary(v);
                std::vector<Cell> flout;
                for (size_t j = 0; j < flin.size(); ++j)
                    flout.push_back(index.at(flin.at(j)));

                faceLists.push_back(flout);
            }

            simpleChainComplex_ = SimpleComplex(dims, scalars, faceLists);
            simpleChainComplexComputed_ = true;
        }
        return simpleChainComplex_;
    }


public:
    MorseData();

    MorseData(ImageData<Value> imageData, float const threshold = -1)
        : imageData_(imageData),
          threshold_(threshold),
          field_(imageData_.complex_),
          fieldComputed_(false),
          criticalComputed_(false),
          criticalSortedComputed_(false),
          basinsComputed_(false),
          watershedsComputed_(false),
          skeletonComputed_(false),
          pathsComputed_(false),
          chainsComputed_(false),
          cochainsComputed_(false),
          simpleChainComplexComputed_(false),
          birthDeathPairsComputed_(false)
    {
    }

    MorseData(ImageData<Value> imageData, std::string const filename)
        : imageData_(imageData),
          threshold_(-1),
          field_(imageData_.xdim(), imageData_.ydim(), imageData_.zdim(),
                 readVolumeData<FieldItem>(filename)),
          fieldComputed_(true),
          criticalComputed_(false),
          criticalSortedComputed_(false),
          basinsComputed_(false),
          watershedsComputed_(false),
          skeletonComputed_(false),
          pathsComputed_(false),
          chainsComputed_(false),
          cochainsComputed_(false),
          simpleChainComplexComputed_(false),
          birthDeathPairsComputed_(false)
    {
    }

    unsigned char cellDirection(Cell const v) const
    {
        return field().V().defined(v) ? field().getDirection(v) : 0;
    }

    Cell associatedMinimum(Cell const v) const
    {
        return basins().get(v);
    }

    bool isOnWatershed(Cell const v) const
    {
        return watersheds().get(v);
    }

    float skeletonValue(Cell const v) const
    {
        return skeleton().get(v);
    }

    bool isOnPath(Cell const v) const
    {
        return paths().get(v);
    }

    size_t chainComplexSize() const
    {
        return critical().size();
    }

    Cell chainComplexCell(size_t const i) const
    {
        return critical().at(i);
    }

    std::vector<Cell> chainComplexBoundary(Cell const v) const
    {
        std::vector<std::pair<Cell, int> > tmp = chains().at(v);
        std::vector<Cell> result;

        for (size_t i = 0; i < tmp.size(); ++i)
        {
            std::pair<Cell, int> p = tmp.at(i);
            for (size_t j = 0; j < p.second; ++j)
                result.push_back(p.first);
        }

        return result;
    }

    std::vector<Cell> chainComplexCoboundary(Cell const v) const
    {
        std::vector<std::pair<Cell, int> > tmp = cochains().at(v);
        std::vector<Cell> result;

        for (size_t i = 0; i < tmp.size(); ++i)
        {
            std::pair<Cell, int> p = tmp.at(i);
            for (size_t j = 0; j < p.second; ++j)
                result.push_back(p.first);
        }

        return result;
    }

    std::vector<std::pair<Cell, Cell> > const& birthDeathPairs() const
    {
        if (not birthDeathPairsComputed_)
        {
            std::vector<Cell> const& sources = criticalSorted();
            size_t const n = sources.size();

            std::vector<Pairing<Cell> > const pairs =
                persistencePairing(simpleChainComplex());

            birthDeathPairs_.reserve(pairs.size());
            for (size_t i = 0; i < pairs.size(); ++i)
            {
                size_t const j = pairs.at(i).partner;

                if (j > i)
                    birthDeathPairs_.push_back(
                        std::pair<Cell, Cell>(
                            sources.at(i),
                            j >= n ? sources.at(i) : sources.at(j)));
            }
        }

        return birthDeathPairs_;
    }

    template<typename T>
    std::vector<std::pair<Cell, T> > weights(std::vector<T> const& w) const
    {
        CubicalComplex const& complex = imageData_.complex_;
        std::vector<Cell> const& sources = criticalSorted();
        size_t const n = sources.size();

        SimpleComplex sc = simpleChainComplex();
        std::vector<Pairing<Cell> > const pairs = persistencePairing(sc);

        Partition<size_t> P(pairs.size());

        std::vector<T> tmp = w;
        std::vector<std::pair<Cell, T> > result;
        for (size_t j = 0; j < n; ++j)
        {
            size_t i = pairs.at(j).partner;
            if (i > j)
                continue;
            T val;

            if (pairs.at(i).dimension == 0)
            {
                std::vector<size_t> faces = sc.cellFaces(j);
                size_t k = faces.at(0) == i ? faces.at(1) : faces.at(0);
                val = tmp.at(P.find(i));
                T added = tmp.at(P.find(k));
                P.unite(i, k);
                tmp.at(P.find(i)) = val + added;
            }
            else
            {
                val = tmp.at(j);
                std::vector<Cell> const& chain = pairs.at(i).chain;
                for (size_t k = 0; k < chain.size(); ++k)
                    val = val + tmp.at(chain.at(k));
                tmp.at(j) = val;
            }
            result.push_back(std::make_pair(sources.at(i), val));
        }

        return result;
    }

    std::vector<std::pair<Cell, int> > const weights() const
    {
        Field::Vectors V = field().V();
        Field::Vectors coV = field().coV();
        CubicalComplex const& complex = imageData_.complex_;
        Facets I(complex.xdim(), complex.ydim(), complex.zdim(), false);
        Facets coI(complex.xdim(), complex.ydim(), complex.zdim(), true);

        std::vector<Set<Cell> > w;
        std::vector<Set<Cell> > unstable;
        std::vector<Cell> const& sources = criticalSorted();
        for (size_t i = 0; i < sources.size(); ++i)
        {
            Cell const cell = sources.at(i);
            if (complex.cellDimension(cell) == 0)
                unstable.push_back(Set<Cell>(unstableSet(cell, coV, coI)));
            else
                unstable.push_back(Set<Cell>(unstableSet(cell, V, I)));
            w.push_back(Set<Cell>(i));
        }

        std::vector<std::pair<Cell, Set<Cell> > > tmp = weights(w);
        std::vector<std::pair<Cell, int> > result;
        for (size_t i = 0; i < tmp.size(); ++i) {
            Cell const cell = tmp.at(i).first;
            std::vector<Cell> const& chains = tmp.at(i).second.elements();

            Set<Cell> aggregate;
            for (size_t j = 0; j < chains.size(); ++j)
                aggregate = aggregate + unstable.at(chains.at(j));

            result.push_back(std::make_pair(cell, aggregate.size()));
        }
        return result;
    }
};



} // namespace diamorse
} // namespace anu_am

#endif // !ANU_AM_DIAMORSE_IMAGEANALYSIS_HPP
