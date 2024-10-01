// automatically generated by the FlatBuffers compiler, do not modify
// swiftlint:disable all
// swiftformat:disable all

import FlatBuffers

public struct flatbuffers_goldens_Galaxy: FlatBufferObject, Verifiable {

  static func validateVersion() { FlatBuffersVersion_23_5_26() }
  public var __buffer: ByteBuffer! { return _accessor.bb }
  private var _accessor: Table

  private init(_ t: Table) { _accessor = t }
  public init(_ bb: ByteBuffer, o: Int32) { _accessor = Table(bb: bb, position: o) }

  private enum VTOFFSET: VOffset {
    case numStars = 4
    var v: Int32 { Int32(self.rawValue) }
    var p: VOffset { self.rawValue }
  }

  public var numStars: Int64 { let o = _accessor.offset(VTOFFSET.numStars.v); return o == 0 ? 0 : _accessor.readBuffer(of: Int64.self, at: o) }
  public static func startGalaxy(_ fbb: inout FlatBufferBuilder) -> UOffset { fbb.startTable(with: 1) }
  public static func add(numStars: Int64, _ fbb: inout FlatBufferBuilder) { fbb.add(element: numStars, def: 0, at: VTOFFSET.numStars.p) }
  public static func endGalaxy(_ fbb: inout FlatBufferBuilder, start: UOffset) -> Offset { let end = Offset(offset: fbb.endTable(at: start)); return end }
  public static func createGalaxy(
    _ fbb: inout FlatBufferBuilder,
    numStars: Int64 = 0
  ) -> Offset {
    let __start = flatbuffers_goldens_Galaxy.startGalaxy(&fbb)
    flatbuffers_goldens_Galaxy.add(numStars: numStars, &fbb)
    return flatbuffers_goldens_Galaxy.endGalaxy(&fbb, start: __start)
  }

  public static func verify<T>(_ verifier: inout Verifier, at position: Int, of type: T.Type) throws where T: Verifiable {
    var _v = try verifier.visitTable(at: position)
    try _v.visit(field: VTOFFSET.numStars.p, fieldName: "numStars", required: false, type: Int64.self)
    _v.finish()
  }
}

public struct flatbuffers_goldens_Universe: FlatBufferObject, Verifiable {

  static func validateVersion() { FlatBuffersVersion_23_5_26() }
  public var __buffer: ByteBuffer! { return _accessor.bb }
  private var _accessor: Table

  private init(_ t: Table) { _accessor = t }
  public init(_ bb: ByteBuffer, o: Int32) { _accessor = Table(bb: bb, position: o) }

  private enum VTOFFSET: VOffset {
    case age = 4
    case galaxies = 6
    var v: Int32 { Int32(self.rawValue) }
    var p: VOffset { self.rawValue }
  }

  public var age: Double { let o = _accessor.offset(VTOFFSET.age.v); return o == 0 ? 0.0 : _accessor.readBuffer(of: Double.self, at: o) }
  public var hasGalaxies: Bool { let o = _accessor.offset(VTOFFSET.galaxies.v); return o == 0 ? false : true }
  public var galaxiesCount: Int32 { let o = _accessor.offset(VTOFFSET.galaxies.v); return o == 0 ? 0 : _accessor.vector(count: o) }
  public func galaxies(at index: Int32) -> flatbuffers_goldens_Galaxy? { let o = _accessor.offset(VTOFFSET.galaxies.v); return o == 0 ? nil : flatbuffers_goldens_Galaxy(_accessor.bb, o: _accessor.indirect(_accessor.vector(at: o) + index * 4)) }
  public static func startUniverse(_ fbb: inout FlatBufferBuilder) -> UOffset { fbb.startTable(with: 2) }
  public static func add(age: Double, _ fbb: inout FlatBufferBuilder) { fbb.add(element: age, def: 0.0, at: VTOFFSET.age.p) }
  public static func addVectorOf(galaxies: Offset, _ fbb: inout FlatBufferBuilder) { fbb.add(offset: galaxies, at: VTOFFSET.galaxies.p) }
  public static func endUniverse(_ fbb: inout FlatBufferBuilder, start: UOffset) -> Offset { let end = Offset(offset: fbb.endTable(at: start)); return end }
  public static func createUniverse(
    _ fbb: inout FlatBufferBuilder,
    age: Double = 0.0,
    galaxiesVectorOffset galaxies: Offset = Offset()
  ) -> Offset {
    let __start = flatbuffers_goldens_Universe.startUniverse(&fbb)
    flatbuffers_goldens_Universe.add(age: age, &fbb)
    flatbuffers_goldens_Universe.addVectorOf(galaxies: galaxies, &fbb)
    return flatbuffers_goldens_Universe.endUniverse(&fbb, start: __start)
  }

  public static func verify<T>(_ verifier: inout Verifier, at position: Int, of type: T.Type) throws where T: Verifiable {
    var _v = try verifier.visitTable(at: position)
    try _v.visit(field: VTOFFSET.age.p, fieldName: "age", required: false, type: Double.self)
    try _v.visit(field: VTOFFSET.galaxies.p, fieldName: "galaxies", required: false, type: ForwardOffset<Vector<ForwardOffset<flatbuffers_goldens_Galaxy>, flatbuffers_goldens_Galaxy>>.self)
    _v.finish()
  }
}
