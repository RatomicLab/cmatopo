-- GetTopoGeom
-- This is a helper to get an object of type TopoGeometry that already exists
-- in the relation table. We need this workaround as it seems there's no way in
-- PostGIS to get a TopoGeometry instance of an already existing TopoGeom, we
-- can only create a new one with the same elements.
CREATE OR REPLACE FUNCTION topology.GetTopoGeom(
    toponame varchar, tg_type integer, layer_id integer, tg_id integer)
RETURNS topology.TopoGeometry AS
$$
DECLARE
  ret topology.TopoGeometry;
BEGIN
  ret.id = tg_id;
  SELECT id FROM topology.topology INTO ret.topology_id WHERE name = toponame;
  ret.layer_id = layer_id;
  ret.type = tg_type;
  RETURN ret;
END
$$
LANGUAGE 'plpgsql' VOLATILE STRICT;
