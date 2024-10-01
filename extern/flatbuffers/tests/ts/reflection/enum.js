// automatically generated by the FlatBuffers compiler, do not modify
/* eslint-disable @typescript-eslint/no-unused-vars, @typescript-eslint/no-explicit-any, @typescript-eslint/no-non-null-assertion */
import * as flatbuffers from 'flatbuffers';
import { EnumVal } from '../reflection/enum-val.js';
import { KeyValue } from '../reflection/key-value.js';
import { Type } from '../reflection/type.js';
export class Enum {
    constructor() {
        this.bb = null;
        this.bb_pos = 0;
    }
    __init(i, bb) {
        this.bb_pos = i;
        this.bb = bb;
        return this;
    }
    static getRootAsEnum(bb, obj) {
        return (obj || new Enum()).__init(bb.readInt32(bb.position()) + bb.position(), bb);
    }
    static getSizePrefixedRootAsEnum(bb, obj) {
        bb.setPosition(bb.position() + flatbuffers.SIZE_PREFIX_LENGTH);
        return (obj || new Enum()).__init(bb.readInt32(bb.position()) + bb.position(), bb);
    }
    name(optionalEncoding) {
        const offset = this.bb.__offset(this.bb_pos, 4);
        return offset ? this.bb.__string(this.bb_pos + offset, optionalEncoding) : null;
    }
    values(index, obj) {
        const offset = this.bb.__offset(this.bb_pos, 6);
        return offset ? (obj || new EnumVal()).__init(this.bb.__indirect(this.bb.__vector(this.bb_pos + offset) + index * 4), this.bb) : null;
    }
    valuesLength() {
        const offset = this.bb.__offset(this.bb_pos, 6);
        return offset ? this.bb.__vector_len(this.bb_pos + offset) : 0;
    }
    isUnion() {
        const offset = this.bb.__offset(this.bb_pos, 8);
        return offset ? !!this.bb.readInt8(this.bb_pos + offset) : false;
    }
    mutate_is_union(value) {
        const offset = this.bb.__offset(this.bb_pos, 8);
        if (offset === 0) {
            return false;
        }
        this.bb.writeInt8(this.bb_pos + offset, +value);
        return true;
    }
    underlyingType(obj) {
        const offset = this.bb.__offset(this.bb_pos, 10);
        return offset ? (obj || new Type()).__init(this.bb.__indirect(this.bb_pos + offset), this.bb) : null;
    }
    attributes(index, obj) {
        const offset = this.bb.__offset(this.bb_pos, 12);
        return offset ? (obj || new KeyValue()).__init(this.bb.__indirect(this.bb.__vector(this.bb_pos + offset) + index * 4), this.bb) : null;
    }
    attributesLength() {
        const offset = this.bb.__offset(this.bb_pos, 12);
        return offset ? this.bb.__vector_len(this.bb_pos + offset) : 0;
    }
    documentation(index, optionalEncoding) {
        const offset = this.bb.__offset(this.bb_pos, 14);
        return offset ? this.bb.__string(this.bb.__vector(this.bb_pos + offset) + index * 4, optionalEncoding) : null;
    }
    documentationLength() {
        const offset = this.bb.__offset(this.bb_pos, 14);
        return offset ? this.bb.__vector_len(this.bb_pos + offset) : 0;
    }
    declarationFile(optionalEncoding) {
        const offset = this.bb.__offset(this.bb_pos, 16);
        return offset ? this.bb.__string(this.bb_pos + offset, optionalEncoding) : null;
    }
    static getFullyQualifiedName() {
        return 'reflection.Enum';
    }
    static startEnum(builder) {
        builder.startObject(7);
    }
    static addName(builder, nameOffset) {
        builder.addFieldOffset(0, nameOffset, 0);
    }
    static addValues(builder, valuesOffset) {
        builder.addFieldOffset(1, valuesOffset, 0);
    }
    static createValuesVector(builder, data) {
        builder.startVector(4, data.length, 4);
        for (let i = data.length - 1; i >= 0; i--) {
            builder.addOffset(data[i]);
        }
        return builder.endVector();
    }
    static startValuesVector(builder, numElems) {
        builder.startVector(4, numElems, 4);
    }
    static addIsUnion(builder, isUnion) {
        builder.addFieldInt8(2, +isUnion, +false);
    }
    static addUnderlyingType(builder, underlyingTypeOffset) {
        builder.addFieldOffset(3, underlyingTypeOffset, 0);
    }
    static addAttributes(builder, attributesOffset) {
        builder.addFieldOffset(4, attributesOffset, 0);
    }
    static createAttributesVector(builder, data) {
        builder.startVector(4, data.length, 4);
        for (let i = data.length - 1; i >= 0; i--) {
            builder.addOffset(data[i]);
        }
        return builder.endVector();
    }
    static startAttributesVector(builder, numElems) {
        builder.startVector(4, numElems, 4);
    }
    static addDocumentation(builder, documentationOffset) {
        builder.addFieldOffset(5, documentationOffset, 0);
    }
    static createDocumentationVector(builder, data) {
        builder.startVector(4, data.length, 4);
        for (let i = data.length - 1; i >= 0; i--) {
            builder.addOffset(data[i]);
        }
        return builder.endVector();
    }
    static startDocumentationVector(builder, numElems) {
        builder.startVector(4, numElems, 4);
    }
    static addDeclarationFile(builder, declarationFileOffset) {
        builder.addFieldOffset(6, declarationFileOffset, 0);
    }
    static endEnum(builder) {
        const offset = builder.endObject();
        builder.requiredField(offset, 4); // name
        builder.requiredField(offset, 6); // values
        builder.requiredField(offset, 10); // underlying_type
        return offset;
    }
    unpack() {
        return new EnumT(this.name(), this.bb.createObjList(this.values.bind(this), this.valuesLength()), this.isUnion(), (this.underlyingType() !== null ? this.underlyingType().unpack() : null), this.bb.createObjList(this.attributes.bind(this), this.attributesLength()), this.bb.createScalarList(this.documentation.bind(this), this.documentationLength()), this.declarationFile());
    }
    unpackTo(_o) {
        _o.name = this.name();
        _o.values = this.bb.createObjList(this.values.bind(this), this.valuesLength());
        _o.isUnion = this.isUnion();
        _o.underlyingType = (this.underlyingType() !== null ? this.underlyingType().unpack() : null);
        _o.attributes = this.bb.createObjList(this.attributes.bind(this), this.attributesLength());
        _o.documentation = this.bb.createScalarList(this.documentation.bind(this), this.documentationLength());
        _o.declarationFile = this.declarationFile();
    }
}
export class EnumT {
    constructor(name = null, values = [], isUnion = false, underlyingType = null, attributes = [], documentation = [], declarationFile = null) {
        this.name = name;
        this.values = values;
        this.isUnion = isUnion;
        this.underlyingType = underlyingType;
        this.attributes = attributes;
        this.documentation = documentation;
        this.declarationFile = declarationFile;
    }
    pack(builder) {
        const name = (this.name !== null ? builder.createString(this.name) : 0);
        const values = Enum.createValuesVector(builder, builder.createObjectOffsetList(this.values));
        const underlyingType = (this.underlyingType !== null ? this.underlyingType.pack(builder) : 0);
        const attributes = Enum.createAttributesVector(builder, builder.createObjectOffsetList(this.attributes));
        const documentation = Enum.createDocumentationVector(builder, builder.createObjectOffsetList(this.documentation));
        const declarationFile = (this.declarationFile !== null ? builder.createString(this.declarationFile) : 0);
        Enum.startEnum(builder);
        Enum.addName(builder, name);
        Enum.addValues(builder, values);
        Enum.addIsUnion(builder, this.isUnion);
        Enum.addUnderlyingType(builder, underlyingType);
        Enum.addAttributes(builder, attributes);
        Enum.addDocumentation(builder, documentation);
        Enum.addDeclarationFile(builder, declarationFile);
        return Enum.endEnum(builder);
    }
}
